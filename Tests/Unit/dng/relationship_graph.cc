/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2017 Reed A. Cartwright
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
 *           Reed A. Cartwright <reed@cartwrig.ht>
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

#define BOOST_TEST_MODULE dng::relationship_graph

#include <dng/relationship_graph.h>

#include <vector>
#include <boost/range/adaptor/transformed.hpp>
#include <functional>

#include "../testing.h"

namespace dng {
struct unittest_dng_relationship_graph {
    GETTERS_FOR_MEMBER_VARIABLES(RelationshipGraph,
        (first_founder)(first_nonfounder)
        (first_somatic)(first_library)(roots)
    )
};

template<typename CharType, typename CharTrait>
inline
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const dng::RelationshipGraph::TransitionType& m) {
    o << (int)m;
    return o;
}


}

using u = dng::unittest_dng_relationship_graph;

using namespace dng;
using namespace dng::detail;
using Sex = dng::Pedigree::Sex;
using dng::detail::make_test_range;

BOOST_AUTO_TEST_CASE(test_inheritance_model_from_string) {
    BOOST_CHECK(inheritance_model("autosomal") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("maternal") == InheritanceModel::Maternal);
    BOOST_CHECK(inheritance_model("paternal") == InheritanceModel::Paternal);
    BOOST_CHECK(inheritance_model("x-linked") == InheritanceModel::XLinked);
    BOOST_CHECK(inheritance_model("y-linked") == InheritanceModel::YLinked);
    BOOST_CHECK(inheritance_model("w-linked") == InheritanceModel::WLinked);
    BOOST_CHECK(inheritance_model("z-linked") == InheritanceModel::ZLinked);
    
    BOOST_CHECK(inheritance_model("mitochondrial") == InheritanceModel::Maternal);
    BOOST_CHECK(inheritance_model("xlinked") == InheritanceModel::XLinked);
    BOOST_CHECK(inheritance_model("ylinked") == InheritanceModel::YLinked);
    BOOST_CHECK(inheritance_model("wlinked") == InheritanceModel::WLinked);
    BOOST_CHECK(inheritance_model("zlinked") == InheritanceModel::ZLinked);

    BOOST_CHECK(inheritance_model("a") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("m") == InheritanceModel::Maternal);
    BOOST_CHECK(inheritance_model("p") == InheritanceModel::Paternal);
    BOOST_CHECK(inheritance_model("x") == InheritanceModel::XLinked);
    BOOST_CHECK(inheritance_model("y") == InheritanceModel::YLinked);
    BOOST_CHECK(inheritance_model("w") == InheritanceModel::WLinked);
    BOOST_CHECK(inheritance_model("z") == InheritanceModel::ZLinked);

    BOOST_CHECK(inheritance_model("A") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("M") == InheritanceModel::Maternal);
    BOOST_CHECK(inheritance_model("P") == InheritanceModel::Paternal);
    BOOST_CHECK(inheritance_model("X") == InheritanceModel::XLinked);
    BOOST_CHECK(inheritance_model("Y") == InheritanceModel::YLinked);
    BOOST_CHECK(inheritance_model("W") == InheritanceModel::WLinked);
    BOOST_CHECK(inheritance_model("Z") == InheritanceModel::ZLinked);

    BOOST_CHECK(inheritance_model("Autosomal") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("AUTOSOMAL") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("au") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("aut") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("auto") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("autos") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("autoso") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("autosom") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("autosoma") == InheritanceModel::Autosomal);

    BOOST_CHECK_THROW(inheritance_model(""), std::invalid_argument);
    BOOST_CHECK_THROW(inheritance_model("unknown"), std::invalid_argument);
    BOOST_CHECK_THROW(inheritance_model("aut0somal"), std::invalid_argument);
    BOOST_CHECK_THROW(inheritance_model("b"), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_inheritance_model_to_string) {
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::Autosomal), "AUTOSOMAL");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::Maternal), "MATERNAL");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::Paternal), "PATERNAL");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::XLinked), "X-LINKED");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::YLinked), "Y-LINKED");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::WLinked), "W-LINKED");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::ZLinked), "Z-LINKED");

    BOOST_CHECK_THROW(to_string(InheritanceModel::Unknown), std::invalid_argument);
    BOOST_CHECK_THROW(to_string(static_cast<InheritanceModel>(-1)), std::invalid_argument);
}

using TransitionType = RelationshipGraph::TransitionType;
using transition_t = RelationshipGraph::transition_t;

struct expected_relationship_graph_t {
    int founder;
    int nonfounder;
    int somatic;
    int library;
    int num_nodes;
    int library_nodes_first;
    int library_nodes_second;
    std::vector<int> roots;
    std::vector<std::string> labels;
    std::vector<std::string> library_names;
    std::vector<int> ploidies;
    std::vector<TransitionType> types;
    std::vector<int> parent1;
    std::vector<int> parent2;
    std::vector<float> length1;
    std::vector<float> length2;
};

BOOST_AUTO_TEST_CASE(test_RelationshipGraph_Construct) {
    const double prec = 2*DBL_EPSILON;

    auto test = [prec](const RelationshipGraph &test, const expected_relationship_graph_t& expected) -> void {
        // Check Node Structure
        BOOST_CHECK_EQUAL(u::first_founder(test), expected.founder);
        BOOST_CHECK_EQUAL(u::first_nonfounder(test), expected.nonfounder);
        BOOST_CHECK_EQUAL(u::first_somatic(test), expected.somatic);
        BOOST_CHECK_EQUAL(u::first_library(test), expected.library);

        BOOST_CHECK_EQUAL(test.num_nodes(), expected.num_nodes);
        BOOST_CHECK_EQUAL(test.library_nodes().first, expected.library_nodes_first);
        BOOST_CHECK_EQUAL(test.library_nodes().second, expected.library_nodes_second);

        CHECK_EQUAL_RANGES(u::roots(test), expected.roots);

        // Check Labels and Ploidies
        CHECK_EQUAL_RANGES(test.labels(), expected.labels);
        CHECK_EQUAL_RANGES(test.library_names(), expected.library_names);
        CHECK_EQUAL_RANGES(test.ploidies(), expected.ploidies);

        auto test_transition_types = boost::adaptors::transform(test.transitions(),
            std::mem_fn(&transition_t::type));
        CHECK_EQUAL_RANGES(test_transition_types, expected.types);

        auto test_transition_parent1 = boost::adaptors::transform(test.transitions(),
            std::mem_fn(&RelationshipGraph::transition_t::parent1));
        CHECK_EQUAL_RANGES(test_transition_parent1, expected.parent1);

        auto test_transition_parent2 = boost::adaptors::transform(test.transitions(),
            std::mem_fn(&RelationshipGraph::transition_t::parent2));
        CHECK_EQUAL_RANGES(test_transition_parent2, expected.parent2);

        auto test_transition_length1 = boost::adaptors::transform(test.transitions(),
            std::mem_fn(&RelationshipGraph::transition_t::length1));
        CHECK_EQUAL_RANGES(test_transition_length1, expected.length1);

        auto test_transition_length2 = boost::adaptors::transform(test.transitions(),
            std::mem_fn(&RelationshipGraph::transition_t::length2));
        CHECK_EQUAL_RANGES(test_transition_length2, expected.length2);
    };

    libraries_t libs = {
        {"DadLb", "MomLb","EveLb", "BobLb"},
        {"DadSm", "MomSm", "EveSm", "BobSm"}
    };

    Pedigree quad_ped;
    quad_ped.AddMember("Dad","0","0",Sex::Male,"DadSm");
    quad_ped.AddMember("Mom","0","0",Sex::Female,"MomSm");
    quad_ped.AddMember("Eve","Dad","Mom",Sex::Female,"EveSm");
    quad_ped.AddMember("Bob","Dad","Mom",Sex::Male,"BobSm");

    {
        // Construct Graph
        RelationshipGraph graph;
        graph.Construct(quad_ped, libs, InheritanceModel::Autosomal, 1e-8, 3e-8, 5e-8);

        expected_relationship_graph_t expected = {
            0, 2, 2, 2, 6, 2, 6, // node information
            {0}, // roots
            {"GL/Dad", "GL/Mom", "LB/DadSm/DadLb", "LB/MomSm/MomLb", "LB/EveSm/EveLb", "LB/BobSm/BobLb"},
            {"DadLb", "MomLb", "EveLb", "BobLb"},
            {2,2,2,2,2,2},
            {TransitionType::Founder, TransitionType::Founder,
             TransitionType::Pair, TransitionType::Pair, TransitionType::Trio,TransitionType::Trio},
            {-1,-1,0,1,0,0}, {-1,-1,-1,-1,1,1},
            {0.0,0.0,8e-8,8e-8,9e-8,9e-8}, {0.0,0.0,0.0,0.0,9e-8,9e-8}
        };

        BOOST_TEST_CONTEXT("graph=quad_graph_autosomal") {
            test(graph, expected);
        }
    }

    {
        // Construct Graph
        RelationshipGraph graph;
        graph.Construct(quad_ped, libs, InheritanceModel::XLinked, 1e-8, 3e-8, 5e-8);

        expected_relationship_graph_t expected = {
            0, 2, 2, 2, 6, 2, 6, // node information
            {0}, // roots
            {"GL/Dad", "GL/Mom", "LB/DadSm/DadLb", "LB/MomSm/MomLb", "LB/EveSm/EveLb", "LB/BobSm/BobLb"},
            {"DadLb", "MomLb", "EveLb", "BobLb"},
            {1,2,1,2,2,1},
            {TransitionType::Founder, TransitionType::Founder,
             TransitionType::Pair, TransitionType::Pair, TransitionType::Trio,TransitionType::Pair},
            {-1,-1,0,1,0,1}, {-1,-1,-1,-1,1,-1},
            {0.0,0.0,8e-8,8e-8,9e-8,9e-8}, {0.0,0.0,0.0,0.0,9e-8,0.0}
        };

        BOOST_TEST_CONTEXT("graph=trio_graph_xlinked") {
            test(graph, expected);
        }
    }

    {
        // Construct Graph
        RelationshipGraph graph;
        graph.Construct(quad_ped, libs, InheritanceModel::YLinked, 1e-8, 3e-8, 5e-8);

        expected_relationship_graph_t expected = {
            0, 1, 1, 1, 3, 1, 3, // node information
            {0}, // roots
            {"GL/Dad", "LB/DadSm/DadLb", "LB/BobSm/BobLb"},
            {"DadLb", "BobLb"},
            {1,1,1},
            {TransitionType::Founder, TransitionType::Pair, TransitionType::Pair},
            {-1,0,0}, {-1,-1,-1},
            {0.0,8e-8,9e-8}, {0.0,0.0,0.0}
        };

        BOOST_TEST_CONTEXT("graph=quad_graph_ylinked") {
            test(graph, expected);
        }
    }
}