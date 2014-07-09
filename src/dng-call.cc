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

#include <vector>
#include <boost/range.hpp>
#include <boost/range/algorithm/generate.hpp>
#include <boost/range/irange.hpp>


enum nodes_e {
	a1, a2, a3, a4,
	b1, b2,
	c1, c2,
	a5, c3,
	d1,
	N
};
const char* name[] = {
	"a1", "a2", "a3", "a4",
	"b1", "b2", "c1"
};
	

typedef std::pair<int, int> Edge;
Edge used_by[] = {
	{a1, a2},
	{a1, b1},
	{a2, b1},
	{a3, a4},
	{a3, b2},
	{a4, b2},
	{b1, b2},
	{b1, c1},
	{b2, c1},
	{b1, c2},
	{b2, c2},
	{a1, a5},
	{a1, c3},
	{a5, c3}
};

struct Family {
	int child;
	int dad;
	int mom;
};

Family ped[] {
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 2, 0, 0 },
	{ 3, 0, 0 },
	{ 4, 0, 0 },
	{ 5, 1, 2 },
	{ 6, 3, 4 },
	{ 7, 5, 6 },
	{ 8, 5, 6 },
	{ 9, 0, 0 },
	{ 10, 1, 9 },
	{ 11, 0, 10 },
	{ 12, 0, 11 },
	{ 13, 0, 11 }
};

struct node {
	int id;
	std::vector<int> families;
};

struct family {
	std::vector<int> members;	
};

std::vector<node> node_store;
std::vector<family> family_store;

namespace boost {
  enum edge_family_t { edge_family };
  enum edge_type_t { edge_type };
  enum vertex_art_t { vertex_art };
  
  BOOST_INSTALL_PROPERTY(edge, family);
  BOOST_INSTALL_PROPERTY(edge, type);
  BOOST_INSTALL_PROPERTY(vertex, art);
}


int main(int argc, char* argv[]) {
	using namespace boost;
	using namespace std;

	typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_art_t, bool>,
		property<edge_family_t, std::size_t,property<edge_type_t, std::size_t>>> graph_t;
	typedef graph_traits < graph_t >::vertex_descriptor vertex_t;
	
	// construct graph from pedigree information
	graph_t g(N);
	for(auto a : ped) {
		// check to see if mom and dad have been seen before
		// TODO: check for parent-child inbreeding
		auto id = edge(a.dad, a.mom, g);
		if(!id.second) {
			id = add_edge(a.dad, a.mom, g);
			put(edge_type, g, id.first, 0);
		}
		// add the meiotic edges
		id = add_edge(a.mom, a.child, g);
		put(edge_type, g, id.first, 1);
		id = add_edge(a.dad, a.child, g);
		put(edge_type, g, id.first, 1);
	}
	clear_vertex(0, g);
		
	auto component = get(edge_family, g);
	auto edtype = get(edge_type, g);
	
	// Calculate the biconnected components and articulation points.
	// Store articulation point status in the graph.
	std::vector<vertex_t> art_points;
	auto ret = biconnected_components(g, component, back_inserter(art_points));
	for(auto a : art_points)
		put(vertex_art, g, a, true);
	
	typedef std::vector<std::vector<graph_traits<graph_t>::edge_descriptor>> 
		component_groups_t;
	component_groups_t groups(ret.first);
  	graph_traits<graph_t>::edge_iterator ei, ei_end;
  	for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
  		groups[component[*ei]].push_back(*ei);
  	}
  	graph_traits<graph_t>::vertex_iterator vi, vi_end;
  	for(boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
  		cout << (char)(*vi + 'A') << " " << get(vertex_art, g, *vi) << "\n";
  	}
  	
  	// Identitify the pivot for each group.
  	// The pivot will be the last art. point that has an edge in
  	// the group.  The pivot of the last group doesn't matter.
  	std::vector<vertex_t> pivot(groups.size());
  	for(auto a : art_points) {
  		graph_traits<graph_t>::out_edge_iterator ei, ei_end;
		for(boost::tie(ei, ei_end) = out_edges(a, g); ei != ei_end; ++ei) {
			// Just overwrite existing value so that the last one wins.
			pivot[component(*ei)] = a;
		}
  	}
  	
  	// Detect Family Structure and pivot positions
  	for(int k = 0; k < groups.size(); ++k) {
  		auto &a = groups[k];
  		// Sort edges based on type and target
  		boost::sort(a, [&](auto x, auto y) -> bool { return
  			(edtype(x) < edtype(y)) &&
  			(target(x, g) < target(y, g)); });
  		// Find the range of the parent types
  		auto pos = boost::find_if(a, [&](auto x) -> bool { return edtype(x) != 0; });
  		size_t num_par = distance(a.begin(),pos);
  		
  		// Check to see what type of graph we have
  		if(num_par == 0) {
  			// We have a parent-child single branch
  			if(a.size() != 1)
  				return 1;
  			int parent = source(*pos, g);
  			int child = target(*pos, g);
  		} else if(num_par == 1) {
  			std::vector<graph_traits<graph_t>::vertex_descriptor> vfam;
  			vfam.push_back(source(a.front(), g)); // Dad
  			vfam.push_back(target(a.front(), g)); // Mom
  			while(pos != a.end()) {
  				vfam.push_back(target(*pos, g)); // Child
  				++pos; ++pos; // child edges come in pairs
  			}
  			auto pos2 = boost::find(vfam, pivot[k]);
  			size_t p = distance(vfam.begin(),pos2);
  			cout << p;
  			for(auto n : vfam) {
  				cout << " " << (char)(n + 'A');
  			}
  			cout << "\n";
  		} else {
  			return 1;
  		}
  		
   	}
  	
  	// Print groups and pivots
   	for(int k = 0; k < groups.size(); ++k) {
   		auto &a = groups[k];
  		cout << "[";
  		for(auto &b : a) {
  			const char *line = (edtype(b) == 1) ? " -> " : " -- ";
  			cout << " " << (char)(source(b, g) + 'A') << line
  				 << (char)(target(b, g) + 'A');
  		}
  		cout << " ]";
  		cout << " " << (char)(pivot[k] + 'A');	  		
  		cout << "\n";
  	}
	
	return 0;
	
/*	typedef std::list<Vertex> MakeOrder;
	MakeOrder make_order;
	boost::topological_sort(g, std::front_inserter(make_order));

	std::cout << "make ordering: \n";
	for (MakeOrder::iterator i = make_order.begin();
		i != make_order.end(); ++i) {
		std::cout << name[*i] << "\n";
		Graph::in_edge_iterator in_begin, in_end;
		for (boost::tie(in_begin, in_end) = in_edges(*i,g); in_begin != in_end; ++in_begin)
		{   
			std::cout << "    " << name[source(*in_begin,g)] << "\n";
		}
	}
		
	std::cout << std::endl;
*/
	
	size_t sz = boost::size(ped);
	std::vector<int> mat(sz*sz);
	std::vector<int> vsort;
	for(int i=1;i<sz;++i) {
		int m = ped[i].mom;
		int d = ped[i].dad;
		int c = ped[i].child;
		mat[m*sz+m] = 1;
		mat[m*sz+d] = 1;
		mat[m*sz+c] = 1;
		mat[d*sz+m] = 1;
		mat[d*sz+d] = 1;
		mat[d*sz+c] = 1;
		mat[c*sz+m] = 1;
		mat[c*sz+d] = 1;
		mat[c*sz+c] = 1;
		vsort.push_back(c);
	}
	//boost::reverse(vsort);
	for(int i=0;i < mat.size(); ++i) {
		std::cout << mat[i];
		if(i % sz == sz-1)
			std::cout << "\n";
	}
	std::cout << "\n";
	
	for(int g=0;g<vsort.size();++g) {
		// Find first element to sum out
		int j;
		for(j=g;j<vsort.size();++j) {
			// Find all neighbors of v[j];
			int v = vsort[j];
			std::vector<int> vn;
			for(int k=1;k<sz;++k) {
				if(mat[v*sz+k] == 1) {
					vn.push_back(k);
				}
			}
			// Check to see if all neighbors are connected to each other
			for(int x=0;x<vn.size();++x) {
				for(int y=0;y<vn.size();++y) {
					if(mat[vn[x]*sz+vn[y]] != 1)
						goto next_loop;
				}
			}
			std::swap(vsort[j],vsort[g]);
			
			for(int x=0;x<vn.size();++x) {
				mat[vn[x]*sz+v] = 0;
			}
			
			for(int i=0;i < mat.size(); ++i) {
				std::cout << mat[i];
				if(i % sz == sz-1)
					std::cout << "\n";
			}
			std::cout << "\n";
			
			break;
		next_loop:
			continue;
		}
		if(j == vsort.size()) {
			std::cout << "not a zero-loop pedigree" << std::endl;
			break;
		}
		
		
	}
	for(auto v : vsort) {
		std::cout << v << std::endl;
	}
	
	return 0;
}

