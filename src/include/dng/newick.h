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
#ifndef DNG_NEWICK_H
#define DNG_NEWICK_H

#include <boost/graph/adjacency_list.hpp>
#include <string>

namespace dng { namespace newick {

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
	boost::property<boost::vertex_name_t, std::string,
	boost::property<boost::vertex_distance_t, float>>,
	boost::no_property
	> tree_t;

typedef boost::graph_traits<tree_t>::vertex_descriptor node_t;
typedef boost::graph_traits<tree_t>::edge_descriptor branch_t;

}};


#endif // DNG_NEWICK_H

