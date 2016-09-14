/*
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
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
#ifndef UNIT_DNG_RELATIONSHIP_GRAPH_HELPER_H_
#define UNIT_DNG_RELATIONSHIP_GRAPH_HELPER_H_


#include <array>
#include <tuple>
#include <vector>

#include <dng/relationship_graph.h>
#include <dng/graph.h>
#include "../boost_test_helper.h"

#define DNG_GL_PREFIX "GL-"

typedef std::tuple<int, int, graph::EdgeType, float> EdgeInfo;

struct EdgeInfo2{
    int source_vertex;
    int target_vertex;
    graph::EdgeType type;
    float edge_length;
};


std::vector<EdgeInfo> extract_edge_info(Graph &pedigree_graph) {

    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto edge_length = get(boost::edge_length, pedigree_graph);
    auto node_index = get(boost::vertex_index, pedigree_graph);
    std::vector<EdgeInfo> edge_info_vector;
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
        EdgeInfo pi = std::make_tuple(node_index[source(*ei, pedigree_graph)],
                                      node_index[target(*ei, pedigree_graph)],
                                      edge_types[*ei], edge_length[*ei]);
//        std::cout <<  static_cast<std::underlying_type<graph::EdgeType>::type>(edge_types[*ei]) << std::cout;
//        std::cout <<  static_cast<std::underlying_type<graph::EdgeType>::type> (graph::EdgeType::Library) << std::endl;
//        std::cout <<  static_cast<int> (graph::EdgeType::Library) << std::endl;
//        std::cout <<  static_cast<int> (edge_types[*ei]) << "\t" << edge_length[*ei] << std::endl;
        edge_info_vector.push_back(pi);
    }
    sort(edge_info_vector.begin(), edge_info_vector.end());
    return edge_info_vector;
}


void boost_check_equal_edge(EdgeInfo expected, EdgeInfo actual){

    BOOST_CHECK_EQUAL(std::get<0>(expected), std::get<0>(actual));
    BOOST_CHECK_EQUAL(std::get<1>(expected), std::get<1>(actual));
    BOOST_CHECK(std::get<2>(expected) == std::get<2>(actual));
    BOOST_CHECK_EQUAL(std::get<3>(expected), std::get<3>(actual));
    BOOST_TEST_MESSAGE("Edge:" << std::get<0>(expected) << "-" << std::get<1>(expected)) ;
}



#endif /* UNIT_DNG_RELATIONSHIP_GRAPH_HELPER_H_ */
