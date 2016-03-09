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

enum edge_family_t { edge_family };
enum vertex_group_t { vertex_group };

BOOST_INSTALL_PROPERTY(edge, length);
BOOST_INSTALL_PROPERTY(edge, type);
BOOST_INSTALL_PROPERTY(vertex, label);

BOOST_INSTALL_PROPERTY(edge, family);
BOOST_INSTALL_PROPERTY(vertex, group);
}

namespace dng {
namespace graph {
enum struct EdgeType : std::size_t {
    Spousal, Meiotic, Mitotic, Library
};

typedef boost::property<boost::vertex_group_t, std::size_t> VertexGroupProp;
typedef boost::property<boost::vertex_label_t, std::string, VertexGroupProp>
        VertexLabelProp;
typedef boost::property<boost::edge_family_t, std::size_t> EdgeFamilyProp;
typedef boost::property<boost::edge_length_t, float, EdgeFamilyProp>
        EdgeLengthProp;
typedef boost::property<boost::edge_type_t, EdgeType, EdgeLengthProp>
        EdgeTypeProp;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
        VertexLabelProp, EdgeTypeProp> Graph;
} // namespace graph

using graph::Graph;
using graph::EdgeType;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

} // namespace dng


#endif // DNG_GRAPH_H

