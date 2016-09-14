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
#include "fixture_read_trio_from_file.h"
#include "relationship_graph_helper.h"


struct FixturePedigree : public ReadTrioFromFile{


    dng::RelationshipGraph relationship_graph;

    FixturePedigree(std::string s = "FixturePedigree") : ReadTrioFromFile(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);
        relationship_graph.Construct(io_pedigree, rgs, arg.mu, arg.mu_somatic, arg.mu_library);
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
BOOST_FIXTURE_TEST_CASE(test_constructor, FixturePedigree ) {

    BOOST_CHECK_EQUAL(5, relationship_graph.num_nodes() );

    auto workspace = relationship_graph.CreateWorkspace();
    BOOST_CHECK_EQUAL(0, workspace.founder_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.founder_nodes.second);
    BOOST_CHECK_EQUAL(0, workspace.germline_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.germline_nodes.second);
    BOOST_CHECK_EQUAL(2, workspace.somatic_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.somatic_nodes.second);
    BOOST_CHECK_EQUAL(2, workspace.library_nodes.first);
    BOOST_CHECK_EQUAL(5, workspace.library_nodes.second);

    auto labels = relationship_graph.labels();

    const std::vector<std::string> expected_labels = {
        "GL-1", // founder 1
        "GL-2", // founder 2
        "LB-NA12878:Solexa-135852",  // lib 1
        "LB-NA12891:Solexa-135851",  // lib 2
        "LB-NA12892:Solexa-135853"   // lib 3
//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878:Solexa-135852   NA12891:Solexa-135851   NA12892:Solexa-135853
    };
    for (int j = 0; j < 5; ++j) {
        BOOST_CHECK_EQUAL(expected_labels[j], labels[j]);
    }
/* Code relate to labels
//#define DNG_GL_PREFIX "GL-"
//#define DNG_SM_PREFIX "SM-" // define also in newick.cc
//#define DNG_LB_PREFIX "LB-"
//    // Add the labels for the germline nodes
//    labels[0] = DNG_GL_PREFIX "unknown";
//    for (size_t i = 1; i < num_members; ++i) {
//        labels[i] = DNG_GL_PREFIX + pedigree.name(i);
//    }
//    for(auto && a : rgs.libraries()) {
//        vertex_t v = add_vertex(pedigree_graph);
//        labels[v] = DNG_LB_PREFIX + a;
//    }
//
//    if(!labels[u].empty()) {
//        labels_.push_back(labels[u]);
//    } else {
//        labels_.push_back(DNG_SM_PREFIX "unnamed_node_" + util::to_pretty(vid));
//    }
*/
    double expected_som_lib = 2e-8;
    double expected_mu_som_lib = 3e-8;

    auto transitions = relationship_graph.transitions();

    std::vector<RelationshipGraph::transition_t> expected_transitions = { {
            RelationshipGraph::TransitionType::Founder, S_MAX, S_MAX, 0, 0}, {
            RelationshipGraph::TransitionType::Founder, S_MAX, S_MAX, 0, 0}, {
            RelationshipGraph::TransitionType::Germline, 0, 1,
                expected_mu_som_lib, expected_mu_som_lib}, {
            RelationshipGraph::TransitionType::Somatic, 0, S_MAX,
                expected_som_lib, 0}, {
            RelationshipGraph::TransitionType::Somatic, 1, S_MAX,
                expected_som_lib, 0}
    };
    //Transition related code
//    arg.mu, arg.mu_somatic, arg.mu_library);
//    if(edge_types[*ei] == EdgeType::Meiotic) {
//        lengths[*ei] *= mu;
//    } else if(edge_types[*ei] == EdgeType::Mitotic) {
//        lengths[*ei] *= mu_somatic;
//    } else if(edge_types[*ei] == EdgeType::Library) {
//        lengths[*ei] *= mu_library;
//    }
//
//    for(std::size_t i = first_founder_; i < first_nonfounder_; ++i) {
//        transitions_[i] = {TransitionType::Founder, static_cast<size_t>(-1), static_cast<size_t>(-1), 0, 0};
//    }
//    transitions_[child] = {tt, parent, static_cast<size_t>(-1), lengths[*pos], 0};
//    transitions_[child] = { TransitionType::Germline, dad, mom, lengths[*pos], lengths[*(pos + 1)]  };

    for (int k = 0; k < 5; ++k) {
        auto expected = expected_transitions[k];
        auto actual = transitions[k];
        BOOST_CHECK(expected.type == actual.type);
        BOOST_CHECK_EQUAL(expected.parent1, actual.parent1);
        BOOST_CHECK_EQUAL(expected.parent2, actual.parent2);
        BOOST_CHECK_CLOSE(expected.length1, actual.length1, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
        BOOST_CHECK_CLOSE(expected.length2, actual.length2, BOOST_CLOSE_PERCENTAGE_THRESHOLD);

    }
}


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


BOOST_FIXTURE_TEST_CASE(test_parse_io_pedigree, ReadTrioFromFile){

    RelationshipGraph relationship_graph;

    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
    BOOST_CHECK_EQUAL(3, relationship_graph.first_nonfounder_);
    BOOST_CHECK_EQUAL(4, relationship_graph.first_somatic_);
    BOOST_CHECK_EQUAL(0, relationship_graph.first_library_);

    Graph pedigree_graph(relationship_graph.first_somatic_);

    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    BOOST_CHECK_EQUAL(7, num_vertices(pedigree_graph));
    BOOST_CHECK_EQUAL(6, num_edges(pedigree_graph));

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    std::vector<EdgeInfo> expected_vector {
            std::make_tuple(1, 2, graph::EdgeType::Spousal, 0),
            std::make_tuple(1, 3, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(1, 4, graph::EdgeType::Mitotic, 1.0),
            std::make_tuple(2, 3, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(2, 5, graph::EdgeType::Mitotic, 1.0),
            std::make_tuple(3, 6, graph::EdgeType::Mitotic, 1.0)
    };
    BOOST_CHECK_EQUAL(expected_vector.size(), edge_vector.size());
    for (int i = 0; i < expected_vector.size(); ++i) {
        boost_check_equal_edge(expected_vector[i], edge_vector[i]);
    }
}


BOOST_FIXTURE_TEST_CASE(test_add_lib_from_rgs, ReadTrioFromFile) {

    Graph pedigree_graph(io_pedigree.member_count());

    RelationshipGraph relationship_graph;
    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    relationship_graph.AddLibrariesFromReadGroups(pedigree_graph, rgs);

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);

    BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
    BOOST_CHECK_EQUAL(3, relationship_graph.first_nonfounder_);
    BOOST_CHECK_EQUAL(4, relationship_graph.first_somatic_);
    BOOST_CHECK_EQUAL(7, relationship_graph.first_library_);
    BOOST_CHECK_EQUAL(10, relationship_graph.num_nodes_);
    BOOST_CHECK_EQUAL(10, num_vertices(pedigree_graph));
    BOOST_CHECK_EQUAL(9, num_edges(pedigree_graph));

    std::vector<EdgeInfo> expected_edges{
        std::make_tuple(1,2, graph::EdgeType::Spousal, 0.0),
        std::make_tuple(1,3, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(1,4, graph::EdgeType::Mitotic, 1.0),
        std::make_tuple(2,3, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(2,5, graph::EdgeType::Mitotic, 1.0),
        std::make_tuple(3,6, graph::EdgeType::Mitotic, 1.0),
        std::make_tuple(4,8, graph::EdgeType::Library, 1.0),
        std::make_tuple(5,9, graph::EdgeType::Library, 1.0),
        std::make_tuple(6,7, graph::EdgeType::Library, 1.0)
    };

    std::vector<std::string> expected_vertex{
        "GL-unknown",
        "GL-1",
        "GL-2",
        "GL-3",
        "SM-NA12891",
        "SM-NA12892",
        "SM-NA12878",
        "LB-NA12878:Solexa-135852",
        "LB-NA12891:Solexa-135851",
        "LB-NA12892:Solexa-135853"
    };

    for (int v = 0; v < relationship_graph.num_nodes_; ++v) {
        BOOST_CHECK_EQUAL(expected_vertex[v], labels[v]);
    }

    BOOST_CHECK_EQUAL(expected_edges.size(), edge_vector.size());
    for (int i = 0; i < expected_edges.size(); ++i) {
        boost_check_equal_edge(expected_edges[i], edge_vector[i]);
    }
}


BOOST_FIXTURE_TEST_CASE(test_update_edge_lengths, ReadTrioFromFile) {

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

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);

    std::vector<EdgeInfo> expected_edges{
        std::make_tuple(1,2, graph::EdgeType::Spousal, 0.0),
        std::make_tuple(1,3, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(1,4, graph::EdgeType::Mitotic, 0.07),
        std::make_tuple(2,3, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(2,5, graph::EdgeType::Mitotic, 0.07),
        std::make_tuple(3,6, graph::EdgeType::Mitotic, 0.07),
        std::make_tuple(4,8, graph::EdgeType::Library, 0.11),
        std::make_tuple(5,9, graph::EdgeType::Library, 0.11),
        std::make_tuple(6,7, graph::EdgeType::Library, 0.11)
    };

    BOOST_CHECK_EQUAL(expected_edges.size(), edge_vector.size());
    for (int i = 0; i < expected_edges.size(); ++i) {
        boost_check_equal_edge(expected_edges[i], edge_vector[i]);
    }
}


BOOST_FIXTURE_TEST_CASE(test_simplify_pedigree, ReadTrioFromFile) {

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

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);

    BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
    BOOST_CHECK_EQUAL(3, relationship_graph.first_nonfounder_);
    BOOST_CHECK_EQUAL(4, relationship_graph.first_somatic_);
    BOOST_CHECK_EQUAL(7, relationship_graph.first_library_);
    BOOST_CHECK_EQUAL(10, relationship_graph.num_nodes_);
    BOOST_CHECK_EQUAL(10, num_vertices(pedigree_graph));
    BOOST_CHECK_EQUAL(5, num_edges(pedigree_graph));

    std::vector<EdgeInfo> expected_edges{
        std::make_tuple(1,2, graph::EdgeType::Spousal, 0.0),
        std::make_tuple(1,7, graph::EdgeType::Meiotic, 0.23),
        std::make_tuple(1,8, graph::EdgeType::Mitotic, 0.18),
        std::make_tuple(2,7, graph::EdgeType::Meiotic, 0.23),
        std::make_tuple(2,9, graph::EdgeType::Mitotic, 0.18)

    };

    BOOST_CHECK_EQUAL(expected_edges.size(), edge_vector.size());
    for (int i = 0; i < expected_edges.size(); ++i) {
        boost_check_equal_edge(expected_edges[i], edge_vector[i]);
    }
}


BOOST_FIXTURE_TEST_CASE(test_update_labels_node_ids, ReadTrioFromFile) {

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

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);

    BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
    BOOST_CHECK_EQUAL(2, relationship_graph.first_nonfounder_);
    BOOST_CHECK_EQUAL(2, relationship_graph.first_somatic_);
    BOOST_CHECK_EQUAL(2, relationship_graph.first_library_);
    BOOST_CHECK_EQUAL(5, relationship_graph.num_nodes_);
    BOOST_CHECK_EQUAL(10, num_vertices(pedigree_graph));
    BOOST_CHECK_EQUAL(5, num_edges(pedigree_graph));

    std::vector<std::size_t> expected_node_ids = {
            S_MAX,
            0, 1,
            S_MAX,S_MAX,S_MAX,S_MAX,
            2, 3, 4};
    std::vector<std::string> expected_labels {
        "GL-1",
        "GL-2",
        "LB-NA12878:Solexa-135852",
        "LB-NA12891:Solexa-135851",
        "LB-NA12892:Solexa-135853"
    };

    boost_check_equal_vector(expected_labels, relationship_graph.labels_);
    boost_check_equal_vector(expected_node_ids, node_ids);
}


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

