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

#include <vector>
#include <boost/range.hpp>
#include <boost/range/algorithm/generate.hpp>
#include <boost/range/irange.hpp>


enum nodes_e {
	a1, a2, a3, a4,
	b1, b2,
	c1,
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
	{b2, c1}
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
	{ 10, 1, 9 }
};


int main(int argc, char* argv[]) {
/*
	using namespace boost;
	typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;

	Graph g(used_by, used_by + sizeof(used_by) / sizeof(Edge), N);
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	typedef std::list<Vertex> MakeOrder;
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

