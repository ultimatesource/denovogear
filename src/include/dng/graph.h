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
#pragma once
#ifndef DNG_GRAPH_H
#define DNG_GRAPH_H

#include <boost/graph/adjacency_list.hpp>
#include <string>

// Install graph properties for pedigree analysis.
namespace boost {
  enum edge_family_t { edge_family };
  enum vertex_group_t { vertex_group };
  enum edge_type_t { edge_type };
  
  BOOST_INSTALL_PROPERTY(edge, family);
  BOOST_INSTALL_PROPERTY(vertex, group);
  BOOST_INSTALL_PROPERTY(edge, type);
}

namespace dng {

namespace detail {
	using namespace boost;
	typedef adjacency_list<vecS, vecS, bidirectionalS,
		property<vertex_group_t, std::size_t,
		property<vertex_name_t, std::string,
		property<vertex_distance_t, float
		>>>,
		property<edge_family_t, std::size_t,
		property<edge_type_t, std::size_t,
		property<edge_weight_t, float
		>>>
		> graph_t;
}
typedef detail::graph_t Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

};


#endif // DNG_GRAPH_H

