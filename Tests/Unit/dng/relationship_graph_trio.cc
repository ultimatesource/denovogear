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


#define BOOST_TEST_MODULE dng::relationship_graph_trio

#include <dng/relationship_graph.h>

#include "../boost_test_helper.h"
#include "fixture_read_test_from_file.h"
#include "relationship_graph_helper.h"


struct FixturePedigree : public ReadTrioFromFile{

    std::string fixture;
    dng::RelationshipGraph relationship_graph;

    FixturePedigree(std::string s = "FixturePedigree")
            : ReadTrioFromFile(), fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);
        relationship_graph.Construct(io_pedigree, rgs, arg.mu, arg.mu_somatic,
                                     arg.mu_library);
    }

    ~FixturePedigree() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }
};

/*

0     1
|     |
---|---
|  |  |
3  2  4

*/
namespace dng {

BOOST_FIXTURE_TEST_CASE(test_pedigree_inspect, FixturePedigree) {

    BOOST_CHECK_EQUAL(5, relationship_graph.num_nodes());

    std::vector<peel::family_members_t> family =
            relationship_graph.family_members_;
    std::vector<peel::family_members_t> expected_family = {
            {1, 4},
            {0, 1, 2},
            {0, 3}
    };

    for (int f = 0; f < expected_family.size(); ++f) {
        boost_check_equal_vector(expected_family[f], family[f]);
    }

    std::vector<decltype(peel::op::NUM)> ops = relationship_graph.peeling_ops_;
    std::vector<decltype(peel::op::NUM)> expected_ops = {peel::op::UP,
            peel::op::TOFATHER, peel::op::UP};
    boost_check_equal_vector(expected_ops, ops);

    std::vector<decltype(peel::op::NUM)> functions_ops =
            relationship_graph.peeling_functions_ops_;
    std::vector<decltype(peel::op::NUM)> expected_functions_ops = {
            peel::op::UPFAST, peel::op::TOFATHERFAST, peel::op::UP};

    boost_check_equal_vector(expected_functions_ops, functions_ops);
}


// BOOST_FIXTURE_TEST_CASE(test_parse_io_pedigree, ReadTrioFromFile){

//     RelationshipGraph relationship_graph;

//     relationship_graph.SetupFirstNodeIndex(io_pedigree);
//     BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
//     BOOST_CHECK_EQUAL(3, relationship_graph.first_nonfounder_);
//     BOOST_CHECK_EQUAL(4, relationship_graph.first_somatic_);
//     BOOST_CHECK_EQUAL(0, relationship_graph.first_library_);

//     Graph pedigree_graph(relationship_graph.first_somatic_);

//     relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
//     BOOST_CHECK_EQUAL(7, num_vertices(pedigree_graph));
//     BOOST_CHECK_EQUAL(6, num_edges(pedigree_graph));

//     std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
//     std::vector<EdgeInfo> expected_vector {
//             std::make_tuple(1, 2, graph::EdgeType::Spousal, 0),
//             std::make_tuple(1, 3, graph::EdgeType::Meiotic, 1.0),
//             std::make_tuple(1, 4, graph::EdgeType::Mitotic, 1.0),
//             std::make_tuple(2, 3, graph::EdgeType::Meiotic, 1.0),
//             std::make_tuple(2, 5, graph::EdgeType::Mitotic, 1.0),
//             std::make_tuple(3, 6, graph::EdgeType::Mitotic, 1.0)
//     };
//     BOOST_CHECK_EQUAL(expected_vector.size(), edge_vector.size());
//     for (int i = 0; i < expected_vector.size(); ++i) {
//         boost_check_equal_edge(expected_vector[i], edge_vector[i]);
//     }
// }

BOOST_FIXTURE_TEST_CASE(test_create_families_info, ReadTrioFromFile) {

    Graph pedigree_graph(io_pedigree.member_count());

    RelationshipGraph relationship_graph;
    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    relationship_graph.AddLibrariesFromReadGroups(pedigree_graph, rgs);
    double expected_mu = 0.05;
    double expected_mu_somatic = 0.07;
    double expected_mu_library = 0.11;
    relationship_graph.UpdateEdgeLengths(pedigree_graph, expected_mu,
                                         expected_mu_somatic,
                                         expected_mu_library);
    relationship_graph.SimplifyPedigree(pedigree_graph);
    std::vector<size_t> node_ids(relationship_graph.num_nodes_, -1);
    relationship_graph.UpdateLabelsNodeIds(pedigree_graph, rgs, node_ids);

    RelationshipGraph::family_labels_t family_labels;//(num_families);
    std::vector<vertex_t> pivots;//(num_families, dummy_index);
    relationship_graph.CreateFamiliesInfo(pedigree_graph, family_labels, pivots);


    std::vector<std::vector<EdgeInfo>> expected_family_labels {
        {std::make_tuple(2, 9, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(1, 2, graph::EdgeType::Spousal, 0.0),
            std::make_tuple(2, 7, graph::EdgeType::Meiotic, 0.23),
            std::make_tuple(1, 7, graph::EdgeType::Meiotic, 0.23)},
        {std::make_tuple(1, 8, graph::EdgeType::Mitotic, 0.18)},
    };
    std::vector<vertex_t> expected_pivots {2, 1, 0};

    boost_check_equal_vector(expected_pivots, pivots);

    BOOST_CHECK_EQUAL(expected_family_labels.size(), family_labels.size());
    for (int i = 0; i < expected_family_labels.size(); ++i) {
//        boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long> x;
        for (int j = 0; j < expected_family_labels[i].size(); ++j) {
            auto f = family_labels[i][j];
            EdgeInfo actual = std::make_tuple(f.m_source, f.m_target,
                    boost::get(boost::edge_type, pedigree_graph, f),
                    boost::get(boost::edge_length, pedigree_graph, f) );

            boost_check_equal_edge(expected_family_labels[i][j], actual);
        }
    }
}


BOOST_FIXTURE_TEST_CASE(test_create_peeling_ops, ReadTrioFromFile) {

    Graph pedigree_graph(io_pedigree.member_count());

    RelationshipGraph relationship_graph;
    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    relationship_graph.AddLibrariesFromReadGroups(pedigree_graph, rgs);
    double expected_mu = 0.05;
    double expected_mu_somatic = 0.07;
    double expected_mu_library = 0.11;
    relationship_graph.UpdateEdgeLengths(pedigree_graph, expected_mu,
                                         expected_mu_somatic,
                                         expected_mu_library);
    relationship_graph.SimplifyPedigree(pedigree_graph);
    std::vector<size_t> node_ids(relationship_graph.num_nodes_, -1);
    relationship_graph.UpdateLabelsNodeIds(pedigree_graph, rgs, node_ids);

    RelationshipGraph::family_labels_t family_labels;//(num_families);
    std::vector<vertex_t> pivots;//(num_families, dummy_index);
    relationship_graph.CreateFamiliesInfo(pedigree_graph, family_labels, pivots);
    relationship_graph.CreatePeelingOps(pedigree_graph, node_ids, family_labels, pivots);

    std::vector<std::size_t> expected_roots {0};

    std::vector<RelationshipGraph::transition_t> expected_transitions { {
            RelationshipGraph::TransitionType::Founder, S_MAX, S_MAX, 0, 0}, {
            RelationshipGraph::TransitionType::Founder, S_MAX, S_MAX, 0, 0}, {
            RelationshipGraph::TransitionType::Germline, 0, 1, 0.23, 0.23}, {
            RelationshipGraph::TransitionType::Somatic, 0, S_MAX, 0.18, 0}, {
            RelationshipGraph::TransitionType::Somatic, 1, S_MAX, 0.18, 0}
    };

    std::vector<decltype(peel::op::NUM)> expected_peeling_ops = {
            peel::op::UP,
            peel::op::TOFATHER,
            peel::op::UP
    };

    std::vector<peel::family_members_t> expected_family_members {
            {1,4},
            {0,1,2},
            {0,3}
    };

    boost_check_equal_vector(expected_roots, relationship_graph.roots_);

    BOOST_CHECK_EQUAL(expected_transitions.size(), relationship_graph.transitions_.size());
    for (int i = 0; i < expected_transitions.size(); ++i) {
        auto exp = expected_transitions[i];
        auto acutal = relationship_graph.transitions_[i];

        BOOST_CHECK(exp.type == acutal.type);
        BOOST_CHECK_EQUAL(exp.parent1, acutal.parent1);
        BOOST_CHECK_EQUAL(exp.parent2, acutal.parent2);
        BOOST_CHECK_CLOSE(exp.length1, acutal.length1, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
        BOOST_CHECK_CLOSE(exp.length2, acutal.length2, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
    }

    boost_check_equal_vector(expected_peeling_ops, relationship_graph.peeling_ops_);

    BOOST_CHECK_EQUAL(expected_family_members.size(), relationship_graph.family_members_.size());
    for (int i = 0; i < expected_family_members.size(); ++i) {
        auto exp = expected_family_members[i];
        auto actual = relationship_graph.family_members_[i];
        boost_check_equal_vector(exp, actual);
    }
}
} // namespace dng

