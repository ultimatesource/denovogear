/*
 * Copyright (c) 2014 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <utility>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/reverse.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/range/algorithm/find.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

#include <vector>
#include <boost/range.hpp>
#include <boost/range/algorithm/generate.hpp>
#include <boost/range/irange.hpp>

#include <boost/tokenizer.hpp>

// http://www.boost.org/development/requirements.html
// http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml

namespace boost {
  enum edge_family_t { edge_family };
  enum edge_type_t { edge_type };
  enum vertex_art_t { vertex_art };
  
  BOOST_INSTALL_PROPERTY(edge, family);
  BOOST_INSTALL_PROPERTY(edge, type);
  BOOST_INSTALL_PROPERTY(vertex, art);
}

namespace dng {

using boost::multi_index_container;
using boost::multi_index::indexed_by;
using boost::multi_index::random_access;
using boost::multi_index::ordered_unique;
using boost::multi_index::identity;

class Pedigree {
public:
	typedef boost::multi_index_container<std::string,
		indexed_by<
			random_access<>,
			ordered_unique<identity<std::string>>
		>> NameContainer;
	typedef std::vector<std::vector<std::string>> DataTable;
	
	Pedigree() {
	}
	
	// Fetch the name of a member of the pedigree
	const std::string& name(std::size_t id) const {
		return names_[id];
	}
	// How many members, including the dummy, are in the pedigree.
	std::size_t member_count() const {
		return names_.size();
	}
	
	// Given the name of an individual return its id.
	// If the name is not valid, return the id of the 0-th
	// individual.
	std::size_t id(const std::string &name) const {
		auto it = names_.get<1>().find(name);
		if(it == names_.get<1>().end())
			return std::size_t(0);
		return names_.project<0>(it) - names_.begin();
	}
	
	// A reference to the data table
	const DataTable & table() const {
		return table_;
	}
	
	//std::istreambuf_iterator<Char> end_of_stream;
	//std::istreambuf_iterator<Char> stream(o);

	// Parse a string-like object into a pedigree
	// TODO: maybe we can just keep pointers to the deliminators in memory
	// TODO: add comment support
	// TODO: warnings for rows that don't have enough elements?
	// TODO: gender checking
	template<typename Obj>
	bool Parse(Obj &o) {
		using namespace boost;
		using namespace std;
		// Construct the tokenizer
		typedef tokenizer<char_separator<char>,
			typename Obj::const_iterator> tokenizer;
		char_separator<char> sep("\t", "\n", keep_empty_tokens);
		tokenizer tokens(o, sep);
		
		// reset the pedigree string table
		table_.reserve(128);
		table_.resize(1);
		table_.back().clear();
		table_.back().reserve(6);
		
		// Work through tokens and build vectors for each row
		for (auto tok_iter = tokens.begin(); tok_iter != tokens.end();
			++tok_iter) {
			if(*tok_iter == "\n") {
				// Resize so that we have six elements
				table_.back().resize(6,"");
				// Push new row onto table
				table_.emplace_back();
				table_.back().reserve(6);
			} else {
				// Add token to the current row
				table_.back().push_back(*tok_iter);
			}
		}
		table_.back().resize(6,"");
		
		// Go through col 1 and pull out child names
		// Add a dummy 0-th pedigree member to handle
		// unknown individuals.
		names_.clear();
		names_.push_back("");
		for( auto &row : table_)
			names_.push_back(row[1]);
		
		return true;
	}
		
protected:
	NameContainer names_;
	DataTable table_;
};

class PedigreePeeler {
public:
	typedef std::vector<std::vector<std::size_t>> family_members_t;

	bool Construct(const Pedigree& pedigree) {
		using namespace boost;
		using namespace std;
		// Graph to hold pedigree information
		typedef adjacency_list<vecS, vecS, undirectedS,
			property<vertex_art_t, bool>,
			property<edge_family_t, std::size_t,property<edge_type_t, std::size_t>>>
			graph_t;
		typedef graph_traits<graph_t>::vertex_descriptor vertex_t;
		typedef graph_traits<graph_t>::edge_descriptor edge_t;

		graph_t pedigree_graph(pedigree.member_count());
		
		// Go through rows and construct the graph
		for(auto &row : pedigree.table()) {
			// check to see if mom and dad have been seen before
			// TODO: check for parent-child inbreeding
			vertex_t child = pedigree.id(row[1]);
			vertex_t dad = pedigree.id(row[2]);
			vertex_t mom = pedigree.id(row[3]);
			auto id = edge(dad, mom, pedigree_graph);
			if(!id.second) {
				id = add_edge(dad, mom, pedigree_graph);
				put(edge_type, pedigree_graph, id.first, 0);
			}
			// add the meiotic edges
			id = add_edge(mom, child, pedigree_graph);
			put(edge_type, pedigree_graph, id.first, 1);
			id = add_edge(dad, child, pedigree_graph);
			put(edge_type, pedigree_graph, id.first, 1);
		}
		// Remove the dummy individual from the graph
		clear_vertex(pedigree.id(""), pedigree_graph);
		
		auto components = get(edge_family, pedigree_graph);
		auto edge_types = get(edge_type, pedigree_graph);
	
		// Calculate the biconnected components and articulation points.
		vector<vertex_t> articulation_vertices;
		auto result = biconnected_components(pedigree_graph, components,
			back_inserter(articulation_vertices));
		// Store articulation point status in the graph.
		for(auto a : articulation_vertices)
			put(vertex_art, pedigree_graph, a, true);
		
		// Determine which edges belong to which components
		typedef vector<vector<graph_traits<graph_t>::edge_descriptor>> 
			component_groups_t;
		component_groups_t groups(result.first);
	  	graph_traits<graph_t>::edge_iterator ei, ei_end;
	  	for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
	  		groups[components[*ei]].push_back(*ei);
	  	}
	  	graph_traits<graph_t>::vertex_iterator vi, vi_end;
	  	for(boost::tie(vi, vi_end) = vertices(pedigree_graph); vi != vi_end; ++vi) {
	  		cout << (char)(*vi + 'A') << " " << get(vertex_art, pedigree_graph, *vi) << "\n";
	  	}
	  	
	  	// Identitify the pivot for each group.
	  	// The pivot will be the last art. point that has an edge in
	  	// the group.  The pivot of the last group doesn't matter.
	  	vector<vertex_t> pivots(groups.size());
	  	for(auto a : articulation_vertices) {
	  		graph_traits<graph_t>::out_edge_iterator ei, ei_end;
			for(tie(ei, ei_end) = out_edges(a, pedigree_graph); ei != ei_end; ++ei) {
				// Just overwrite existing value so that the last one wins.
				pivots[components[*ei]] = a;
			}
	  	}
	  	// The last pivot is special.
	  	pivots.back() = 0;
	  	
	  	// Reset Family Information
	  	families_.clear();
	  	families_.reserve(128);
	  	
	  	// Detect Family Structure and pivot positions
	  	for(std::size_t k = 0; k < groups.size(); ++k) {
	  		auto &family_edges = groups[k];
	  		// Sort edges based on type and target
	  		boost::sort(family_edges, [&](edge_t x, edge_t y) -> bool { return
	  			(edge_types(x) < edge_types(y)) &&
	  			(target(x, pedigree_graph) < target(y, pedigree_graph)); });
	  		
	  		// Find the range of the parent types
	  		auto pos = boost::find_if(family_edges, [&](edge_t x) -> bool { 
	  			return edge_types(x) != 0; });
	  		size_t num_parent_edges = distance(family_edges.begin(),pos);
	  		
	  		
	  		// Check to see what type of graph we have
	  		if(num_parent_edges == 0) {
	  			// If we do not have a parent-child single branch,
	  			// we can't construct the pedigree.
	  			// TODO: Write Error message
	  			if(family_edges.size() != 1)
	  				return false;
	  			vertex_t parent = source(*pos, pedigree_graph);
	  			vertex_t child = target(*pos, pedigree_graph);
	  			// TODO: What do we do here?
	  		} else if(num_parent_edges == 1) {
	  			// We have a nuclear family with 1 or more children
	  			families_.emplace_back();
	  			families_.back().push_back(source(family_edges.front(), pedigree_graph)); // Dad
	  			families_.back().push_back(target(family_edges.front(), pedigree_graph)); // Mom
	  			while(pos != family_edges.end()) {
	  				families_.back().push_back(target(*pos, pedigree_graph)); // Child
	  				++pos; ++pos; // child edges come in pairs
	  			}
	  			auto pivot_pos = boost::find(families_.back(), pivots[k]);
	  			size_t p = distance(families_.back().begin(),pivot_pos);
	  			if(pivots[k] == 0 ) {
	  				p = 0;
	  			} else if(p > 2) {
	  				swap(families_.back()[p], families_.back()[2]);
	  				p = 2;
	  			}
	  			// TODO: Assert that p makes sense

	  			cout << p;
	  			for(auto n : families_.back()) {
	  				cout << " " << (char)(n + 'A');
	  			}
	  			cout << "\n";
	  		} else {
	  			// Not a zero-loop pedigree
	  			// TODO: write error message
	  			return false;
	  		}
	   	}
		return true;
	}
	
protected:
	family_members_t families_;
};

}; // namespace dng

int main(int argc, char* argv[]) {
	using namespace boost;
	using namespace std;
	
	std::string str =
		"1\t1\t0\t0\n"
		"1\t2\t0\t0\n"
		"1\t3\t1\t2\n"
		"1\t4\t1\t2\n"
		"1\t5\t0\t0\n"
		"1\t6\t5\t2\n"
	;
	dng::Pedigree pedi;
	pedi.Parse(str);
	dng::PedigreePeeler peel;
	peel.Construct(pedi);
		
	
	return 0;
}

