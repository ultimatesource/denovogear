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
	{ 11, 0, 10 }
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

int build_graph(Family *p, size_t sz) {
	std::map<std::pair<int,int>, int> fam_map;
	// reset node_store
	node_store.resize(sz);
	for(size_t k=0;k<sz;++k) {
		int u = p[k].child;
		node_store[u].id = u;
		node_store[u].families.clear();
	}
	// build graph of families
	family_store.clear();
	for(size_t k=0;k<sz;++k) {
		if(p[k].dad == 0 || p[k].mom == 0)
			continue;
		// check to see if parent pair already exists
		auto key = std::make_pair(p[k].dad, p[k].mom);
		auto ret = fam_map.find(key);
		if(ret != fam_map.end()) {
			// family exists, so we just add the child and add it to the family
			family_store[ret->second].members.push_back(p[k].child);
			node_store[p[k].child].families.push_back(family_store.size()-1);
		} else {
			// family doesn't exist, so we construct it
			family_store.push_back(family());
			family_store.back().members.assign({p[k].mom,p[k].dad,p[k].child});
			int fam_id = family_store.size()-1;
			fam_map.emplace(key, fam_id);
			// add the family ids to the nodes
			node_store[p[k].mom].families.push_back(fam_id);
			node_store[p[k].dad].families.push_back(fam_id);
			node_store[p[k].child].families.push_back(fam_id);
		}
	}
	return 0;
}

int construct_peeling_order() {
	size_t sz = family_store.size();
	std::vector<int> peel_order;
	std::vector<int> peeled(sz,0);
	// find all families that are peelable
	for(int i=0;i<sz;++i) {
		auto &fam = family_store[i];
	}
	return 0;
}

namespace boost {
  enum edge_family_t { edge_family };
  enum edge_type_t { edge_type };
  
  BOOST_INSTALL_PROPERTY(edge, family);
  BOOST_INSTALL_PROPERTY(edge, type);
}


int main(int argc, char* argv[]) {
	using namespace boost;
	using namespace std;

	typedef adjacency_list<vecS, vecS, undirectedS,no_property,
		property<edge_family_t, std::size_t,property<edge_type_t, std::size_t>>> graph_t;
	typedef graph_traits < graph_t >::vertex_descriptor vertex_t;
	std::vector<vertex_t> art_points;

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
	auto ret = biconnected_components(g, component,std::back_inserter(art_points));
	
	typedef std::vector<std::vector<graph_traits<graph_t>::edge_iterator>> 
		component_groups_t;
	component_groups_t groups(ret.first);
  	graph_traits<graph_t>::edge_iterator ei, ei_end;
  	for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
  		groups[component[*ei]].push_back(ei);
  	}
  	auto art = art_points.begin();
  	for(auto &a : groups) {
  		cout << "[";
  		for(auto &b : a) {
  			const char *line = (edtype(*b) == 1) ? " -> " : " -- ";
  			cout << " " << (char)(source(*b, g) + 'A') << line
  				 << (char)(target(*b, g) + 'A');
  		}
  		cout << " ]";
  		// check to see if there is only one spousal pair
  		ei = ei_end;
  		for(auto &b : a) {
  			if(edtype(*b) != 0)
  				continue;
  			if(ei != ei_end)
  				return 0;
  			ei = b;
  		}
  		
  		
  		if(art != art_points.end()) {
  			cout << " " << (char)((*art)+'A');
  			++art;
  		}
  		
  		
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

