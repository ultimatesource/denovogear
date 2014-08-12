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
enum edge_length_t { edge_length };
enum edge_type_t { edge_type };
enum vertex_label_t { vertex_label };

BOOST_INSTALL_PROPERTY(edge, length);
BOOST_INSTALL_PROPERTY(edge, type);
BOOST_INSTALL_PROPERTY(vertex, label);
}

namespace dng { namespace graph {
typedef boost::property<boost::vertex_label_t, std::string> VertexLabel;
typedef boost::property<boost::edge_type_t, std::size_t> EdgeType;
typedef boost::property<boost::edge_length_t, float, EdgeType> EdgeLength;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
	VertexLabel, EdgeLength> Graph;
}

using graph::Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

};


#endif // DNG_GRAPH_H

