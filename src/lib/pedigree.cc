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

#include <dng/pedigree.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/range/algorithm/find.hpp>

bool dng::Pedigree::Initialize(double theta, double mu) {
	using namespace Eigen;
	using namespace std;
	// Construct Genotype Prior
	double nuc_freq[4] = { 0.25, 0.25, 0.25, 0.25 };
	double alpha[4] = {
		theta*nuc_freq[0], theta*nuc_freq[1],
		theta*nuc_freq[2], theta*nuc_freq[3]
	};
	// Use a parent-independent mutation model, which produces a
	// beta-binomial
	genotype_prior_ <<
	        alpha[0]*(1.0+alpha[0])/theta/(1.0+theta), // AA
	    2.0*alpha[0]*(    alpha[1])/theta/(1.0+theta), // AC
	    2.0*alpha[0]*(    alpha[2])/theta/(1.0+theta), // AG
	    2.0*alpha[0]*(    alpha[3])/theta/(1.0+theta), // AT
	        alpha[1]*(1.0+alpha[1])/theta/(1.0+theta), // CC
	    2.0*alpha[1]*(    alpha[2])/theta/(1.0+theta), // CG
	    2.0*alpha[1]*(    alpha[3])/theta/(1.0+theta), // CT
	        alpha[2]*(1.0+alpha[2])/theta/(1.0+theta), // GG
	    2.0*alpha[2]*(    alpha[3])/theta/(1.0+theta), // GT
	        alpha[3]*(1.0+alpha[3])/theta/(1.0+theta); // GG
			
	// Construct Mutation Process
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
	for(int i = 0; i < 10; ++i) {
		for(int j = 0; j < 10; ++j) {
			int p = 10*i+j;
			for(int k = 0; k < 10; ++k) {
				meiosis_(p,k) =
					  ( m(nucleotides[i][0],nucleotides[k][0])
					  + m(nucleotides[i][1],nucleotides[k][0]))
				    * ( m(nucleotides[j][0],nucleotides[k][1])
				      + m(nucleotides[j][1],nucleotides[k][1]))/4.0 ;
				if(nucleotides[k][0] == nucleotides[k][1])
					continue;
				meiosis_(p,k) +=
					  ( m(nucleotides[i][0],nucleotides[k][1])
					  + m(nucleotides[i][1],nucleotides[k][1]))
				    * ( m(nucleotides[j][0],nucleotides[k][0])
				      + m(nucleotides[j][1],nucleotides[k][0]))/4.0 ;
			}
		}
	}
			
	return true;
}

bool dng::Pedigree::Construct(const io::Pedigree& pedigree) {
	using namespace boost;
	using namespace std;
	//using dng::newick::vertex_t;
	//using dng::newick::edge_t;
	//using dng::newick::Graph;
	//typedef Graph graph_t;
	
  	// Reset Family Information
  	roots_.clear();
  	family_members_.clear();
  	family_members_.reserve(128);
  	peeling_op_.clear();
  	peeling_op_.reserve(128);
			
	num_members_ = pedigree.member_count();
	
	Graph pedigree_graph(num_members_);
	graph_traits<Graph>::edge_iterator ei, ei_end;
	graph_traits<Graph>::vertex_iterator vi, vi_end;
	
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
	
	auto groups = get(vertex_group, pedigree_graph);
	auto families = get(edge_family, pedigree_graph);
	auto edge_types = get(edge_type, pedigree_graph);
	
	// Calculate the connected components.  This defines independent sections
	// of the graph.
	std::size_t num_groups = connected_components(pedigree_graph, groups);
	
	// Calculate the biconnected components and articulation points.
	// This defines "nuclear" families and pivot indiviuduals.
	// Nodes which have no edges will not be part of any family.
	vector<vertex_t> articulation_vertices;
	std::size_t num_families = biconnected_components(pedigree_graph, families,
		back_inserter(articulation_vertices)).first;
			
	// Determine which edges belong to which nuclear families.
	typedef vector<vector<graph_traits<Graph>::edge_descriptor>> 
		family_labels_t;
	family_labels_t family_labels(num_families);
  	for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei)
  		family_labels[families[*ei]].push_back(*ei);
	
	// Determine the last family in each group.  All singleton groups will have
	// a value of -1 since they have no family assignment.
	typedef deque<std::size_t> root_families_t;
	root_families_t root_families(num_groups,-1);  	
  	for(std::size_t f = 0; f < family_labels.size(); ++f) {
  		// last one wins
  		auto first_edge = family_labels[f][0];
  		auto src_vertex = source(first_edge, pedigree_graph);
		root_families[groups(src_vertex)] = f;
  	}
  	
  	// Identitify the pivot for each family.
  	// The pivot will be the last art. point that has an edge in
  	// the group.  The pivot of the last group doesn't matter.
  	vector<vertex_t> pivots(num_families);
  	for(auto a : articulation_vertices) {
  		graph_traits<Graph>::out_edge_iterator ei, ei_end;
		for(tie(ei, ei_end) = out_edges(a, pedigree_graph); ei != ei_end; ++ei) {
			// Just overwrite existing value so that the last one wins.
			pivots[families[*ei]] = a;
		}
  	}
  	// Root Pivots are special
  	for(auto f : root_families) {
  		if(f == -1) // skip all singleton groups
  			continue;
  		pivots[f] = 0;
  	}

	// Non-dummy singleton groups are root elements.
	for(tie(vi, vi_end) = vertices(pedigree_graph); vi != vi_end; ++vi) {
		if(degree(*vi,pedigree_graph) > 0 || *vi == 0)
			continue;
		roots_.push_back(*vi);
	}

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
  			// A family without a pivot is a root family
  			if(pivots[k] == 0 ) {
  				p = 0;
  				roots_.push_back(family_members_.back()[0]);
  			} else if(p > 2) {
  				swap(family_members_.back()[p], family_members_.back()[2]);
  				p = 2;
  			}
  			// TODO: Assert that p makes sense
			
			switch(p) {
			case 0:
				peeling_op_.emplace_back(&Op::PeelToFather);
				break;
			case 1:
				peeling_op_.emplace_back(&Op::PeelToMother);
				break;
			case 2:
			default:
				peeling_op_.emplace_back(&Op::PeelToChild);
				break;
			};
  		} else {
  			// Not a zero-loop pedigree
  			// TODO: write error message
  			return false;
  		}
   	}
	return true;
}


