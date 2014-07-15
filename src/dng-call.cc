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

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/KroneckerProduct>

// http://www.boost.org/development/requirements.html
// http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml

namespace boost {
  enum edge_family_t { edge_family };
  enum edge_group_t { edge_group };
  enum edge_type_t { edge_type };
  enum vertex_art_t { vertex_art };
  
  BOOST_INSTALL_PROPERTY(edge, family);
  BOOST_INSTALL_PROPERTY(edge, group);
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

	typedef void (PedigreePeeler::*PeelOp)(std::size_t);

	typedef Eigen::Array<double, 10, 1> Vector10d;
	typedef Eigen::Matrix<double, 10, 10> Matrix10d;
	typedef Eigen::Matrix<double, 10, 10, Eigen::RowMajor> RowMatrix10d;
	typedef Eigen::Array<double, 100, 1> Vector100d;
	typedef Eigen::Matrix<double, 100, 100> Matrix100d;
	
	typedef Eigen::Matrix<double, 100, 10> MeiosisMatrix;

	typedef std::vector<Vector10d, Eigen::aligned_allocator<Vector10d>> IndividualBuffer;

	bool Initialize(double mu, double theta) {
		using namespace Eigen;
		using namespace std;
		// Construct Genotype Prior
		double p_hom = 1.0/(1.0+theta)/4.0;
		double p_het = theta/(1.0+theta)/6.0;
		genotype_prior_ <<
			p_hom, p_het, p_het, p_het,
			       p_hom, p_het, p_het,
			              p_hom, p_het,
			                     p_hom
		;
		double nuc_freq[4] = { 0.25, 0.25, 0.25, 0.25 };
		double beta = 1.0;
		for(auto d : nuc_freq)
			beta -= d*d;
		beta = 1.0/beta;
		beta = exp(-beta*mu);
		Eigen::Matrix4d m;
		for(int i : {0,1,2,3}) {
			for(int j : {0,1,2,3}) {
				m(i,j) = nuc_freq[i]*(1.0-beta);
			}
			m(i,i) += beta;
		}
		int nuc1[10] = {0,0,0,0,1,1,1,2,2,3};
		int nuc2[10] = {0,1,2,3,1,2,3,2,3,3};
		for(int i = 0; i < 10; ++i) {
			for(int j = 0; j < 10; ++j) {
				int p = 10*i+j;
				for(int k = 0; k < 10; ++k) {
					meiosis_(p,k) =
						  ( m(nuc1[i],nuc1[k]) + m(nuc2[i],nuc1[k]) )
					    * ( m(nuc1[j],nuc2[k]) + m(nuc2[j],nuc2[k]) )/4.0 ;
					if(nuc1[k] == nuc2[k])
						continue;
					meiosis_(p,k) +=
					      ( m(nuc1[i],nuc2[k]) + m(nuc2[i],nuc2[k]) )
					    * ( m(nuc1[j],nuc1[k]) + m(nuc2[j],nuc1[k]) )/4.0 ;
				}
			}
		}
				
		return true;
	}

	bool Construct(const Pedigree& pedigree) {
		// TODO: Check for independent families
		using namespace boost;
		using namespace std;
		// Graph to hold pedigree information
		typedef adjacency_list<vecS, vecS, undirectedS,
			property<vertex_art_t, bool>,
			property<edge_family_t, std::size_t,
			property<edge_group_t, std::size_t, 
			property<edge_type_t, std::size_t
			>>>
			> graph_t;
		typedef graph_traits<graph_t>::vertex_descriptor vertex_t;
		typedef graph_traits<graph_t>::edge_descriptor edge_t;
		
		num_members_ = pedigree.member_count();
		
		graph_t pedigree_graph(num_members_);
		graph_traits<graph_t>::edge_iterator ei, ei_end;
		
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
		
		//auto groups = get(edge_groups, pedigree_graph);
		auto families = get(edge_family, pedigree_graph);
		auto edge_types = get(edge_type, pedigree_graph);
		
		// Calculate the connected componetns
		std::size_t num_groups = connected_components(pedigree_graph, groups);
		
		// Calculate the biconnected components and articulation points.
		vector<vertex_t> articulation_vertices;
		std::size_t num_families = biconnected_components(pedigree_graph, families,
			back_inserter(articulation_vertices)).first;

		// Store articulation point status in the graph.
		for(auto a : articulation_vertices)
			put(vertex_art, pedigree_graph, a, true);

		// Determine the last family in each group
		typedef vector<std::size_t> last_family_in_group_t;
		last_family_in_group_t last_family_in_group(num_groups,0);  	
	  	for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
	  		last_family_in_group[groups[*ei]] = std::max(families[*ei],
	  			last_family_in_group[groups[*ei]]);
	  	}

		// Determine which edges belong to which nuclear families
		typedef vector<vector<graph_traits<graph_t>::edge_descriptor>> 
			family_labels_t;
		family_labels_t family_labels(num_families);
	  	for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
	  		family_labels[families[*ei]].push_back(*ei);
	  	}
	  	
	  	// Identitify the pivot for each family.
	  	// The pivot will be the last art. point that has an edge in
	  	// the group.  The pivot of the last group doesn't matter.
	  	vector<vertex_t> pivots(num_families);
	  	for(auto a : articulation_vertices) {
	  		graph_traits<graph_t>::out_edge_iterator ei, ei_end;
			for(tie(ei, ei_end) = out_edges(a, pedigree_graph); ei != ei_end; ++ei) {
				// Just overwrite existing value so that the last one wins.
				pivots[families[*ei]] = a;
			}
	  	}
	  	// The last pivot is special.
	  	pivots.back() = 0;
	  	
	  	// Reset Family Information
	  	family_members_.clear();
	  	family_members_.reserve(128);
	  	peeling_op_.clear();
	  	peeling_op_.reserve(128);
	  	
	  	// Detect Family Structure and pivot positions
	  	for(std::size_t k = 0; k < family_labels.size(); ++k) {
	  		auto &family_edges = family_labels[k];
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
	  			family_members_.emplace_back();
	  			family_members_.back().push_back(source(family_edges.front(), pedigree_graph)); // Dad
	  			family_members_.back().push_back(target(family_edges.front(), pedigree_graph)); // Mom
	  			while(pos != family_edges.end()) {
	  				family_members_.back().push_back(target(*pos, pedigree_graph)); // Child
	  				++pos; ++pos; // child edges come in pairs
	  			}
	  			auto pivot_pos = boost::find(family_members_.back(), pivots[k]);
	  			size_t p = distance(family_members_.back().begin(),pivot_pos);
	  			// A family without a pivot is the final family
	  			if(pivots[k] == 0 ) {
	  				p = 0;
	  				root_ = family_members_.back()[0];
	  			} else if(p > 2) {
	  				swap(family_members_.back()[p], family_members_.back()[2]);
	  				p = 2;
	  			}
	  			// TODO: Assert that p makes sense
				
				switch(p) {
				case 0:
					peeling_op_.push_back(&PedigreePeeler::PeelToFather);
					break;
				case 1:
					peeling_op_.push_back(&PedigreePeeler::PeelToMother);
					break;
				case 2:
				default:
					peeling_op_.push_back(&PedigreePeeler::PeelToChild);
					break;
				};
				
	  			cout << p;
	  			for(auto n : family_members_.back()) {
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
	
	double CalculateLogLikelihood(const IndividualBuffer &penetrances) {
		// TODO: Assert that length of penetrances is equal to num_members_;
		// Copy Penetrance values into the lower buffer
		// TODO: Eliminate this copy????
		lower_ = penetrances;
		// Copy genotype Priors
		upper_.assign(num_members_, genotype_prior_);
		
		// Peel Pedigree on family at a time
		for(std::size_t i = 0; i < peeling_op_.size(); ++i)
			(this->*(peeling_op_[i]))(i);
			
		// Sum over root
		return (lower_[root_] * upper_[root_]).sum();
	}
	
protected:
	family_members_t family_members_;
	std::size_t num_members_;
	std::size_t root_;
	
	IndividualBuffer upper_; // Holds P(Data & G=g)
	IndividualBuffer lower_; // Holds P(Data | G=g)
	
	Vector10d genotype_prior_; // Holds P(G | theta)

	MeiosisMatrix meiosis_;

	std::vector<PeelOp> peeling_op_;

	void PeelToFather(std::size_t id) {
		using namespace Eigen;
		auto family_members = family_members_[id];
		Vector100d buffer;
		buffer.setOnes();
		// Sum over children
		for(std::size_t i = 2; i < family_members.size(); i++) {
			buffer *= (meiosis_ * lower_[family_members[i]].matrix()).array();
		}
		// Include Mom
		Map<Matrix10d, Aligned> mat(buffer.data());
		lower_[family_members[0]] *= (mat *
			(upper_[family_members[1]]*lower_[family_members[1]]).matrix()).array();
	}

	void PeelToMother(std::size_t id) {
		using namespace Eigen;
		auto family_members = family_members_[id];
		Vector100d buffer;
		buffer.setOnes();
		// Sum over children
		for(std::size_t i = 2; i < family_members.size(); i++) {
			buffer *= (meiosis_ * lower_[family_members[i]].matrix()).array();
		}
		// Include Dad
		Map<RowMatrix10d, Aligned> mat(buffer.data());
		lower_[family_members[1]] *= (mat *
			(upper_[family_members[0]]*lower_[family_members[0]]).matrix()).array();
	}

	void PeelToChild(std::size_t id) {
		using namespace Eigen;
		auto family_members = family_members_[id];
		Vector100d buffer;
		buffer.setOnes();
		// Sum over children
		for(std::size_t i = 3; i < family_members.size(); i++) {
			buffer *= (meiosis_ * lower_[family_members[i]].matrix()).array();
		}
		// Parents
		buffer *= kroneckerProduct(
			(lower_[family_members[0]] * upper_[family_members[0]]).matrix(),
			(lower_[family_members[1]] * upper_[family_members[1]]).matrix()
		).array();
		
		upper_[family_members[2]].matrix() = meiosis_.transpose() * buffer.matrix();
	}
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
		"1\t6\t4\t5\n"
	;
	dng::Pedigree pedi;
	pedi.Parse(str);
	dng::PedigreePeeler peel;
	peel.Initialize(1e-12,1e-12);
	peel.Construct(pedi);
	
	dng::PedigreePeeler::IndividualBuffer buf;
	buf.resize(7, dng::PedigreePeeler::Vector10d::Zero());
	buf[1][0] = 1.0;
	buf[2][0] = 1.0;
	buf[3][0] = 1.0;
	buf[4][0] = 1.0;
	buf[5][0] = 1.0;
	buf[6][0] = 1.0;
	
	double d = peel.CalculateLogLikelihood(buf);
	cout << d << endl;
	
	return 0;
}

