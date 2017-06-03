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

struct topology_t {
    std::string name;
    libraries_t libs;
    Pedigree ped;
};

struct rates_t {
    double germline;
    double somatic;
    double library;
};

template<typename F>
void run_graph_tests(F test, double prec = 2*DBL_EPSILON) {
    topology_t trio{"trio"};
    trio.libs = {
        {"Mom", "Dad", "Eve"},
        {"Mom", "Dad", "Eve"}
    };
    trio.ped.AddMember("Dad","0","0",Sex::Male,"");
    trio.ped.AddMember("Mom","0","0",Sex::Female,"");
    trio.ped.AddMember("Eve","Dad","Mom",Sex::Female,"");

    test(trio, {1e-8,1e-8,1e-8}, InheritanceModel::Autosomal, prec);
}

BOOST_AUTO_TEST_CASE(test_trio_RelationshipGraph_Construct) {
    const double prec = 2*DBL_EPSILON;

    libraries_t libs = {
        {"MomLb", "DadLb", "EveLb"},
        {"MomSm", "DadSm", "EveSm"}
    };

    Pedigree ped;
    ped.AddMember("Dad","0","0",Sex::Male,"DadSm");
    ped.AddMember("Mom","0","0",Sex::Female,"MomSm");
    ped.AddMember("Eve","Dad","Mom",Sex::Female,"EveSm");

    // Construct Graph
    RelationshipGraph g;
    g.Construct(ped, libs, InheritanceModel::Autosomal, 1e-8, 1e-8, 1e-8);

    // Check Node Structure
    BOOST_CHECK_EQUAL(u::first_founder(g), 0);
    BOOST_CHECK_EQUAL(u::first_nonfounder(g), 2);
    BOOST_CHECK_EQUAL(u::first_somatic(g), 2);
    BOOST_CHECK_EQUAL(u::first_library(g), 2);

    BOOST_CHECK_EQUAL(g.num_nodes(), 5);
    BOOST_CHECK_EQUAL(std::get<0>(g.library_nodes()), 2);
    BOOST_CHECK_EQUAL(std::get<1>(g.library_nodes()), 5);

    std::vector<int> expected_roots = {0};
    CHECK_EQUAL_RANGES(u::roots(g), expected_roots);

    // Check Labels and Ploidies
    std::vector<std::string> expected_labels = { 
        "GL/Dad", "GL/Mom", "LB/MomSm/MomLb", "LB/DadSm/DadLb", "LB/EveSm/EveLb"
    };
    CHECK_EQUAL_RANGES(g.labels(), expected_labels);

    std::vector<std::string> expected_library_names = {
        "MomLb", "DadLb", "EveLb"
    };
    CHECK_EQUAL_RANGES(g.library_names(), expected_library_names);

    std::vector<int> expected_ploidies(5,2);
    CHECK_EQUAL_RANGES(g.ploidies(), expected_ploidies);

    // Check Transitions
    using TransitionType = RelationshipGraph::TransitionType;
    using transition_t = RelationshipGraph::transition_t;

    std::vector<TransitionType> expected_transition_types = {
       TransitionType::Founder, TransitionType::Founder,
       TransitionType::Pair, TransitionType::Pair, TransitionType::Trio
    };
    std::vector<int> expected_transition_parent1 = {-1,-1,1,0,0};
    std::vector<int> expected_transition_parent2 = {-1,-1,-1,-1,1};
    std::vector<float> expected_transition_length1 = {0.0,0.0,2e-8,2e-8,3e-8};
    std::vector<float> expected_transition_length2 = {0.0,0.0,0.0,0.0,3e-8};

    auto test_transition_types = boost::adaptors::transform(g.transitions(),
        std::mem_fn(&transition_t::type));
    CHECK_EQUAL_RANGES(test_transition_types, expected_transition_types);

    auto test_transition_parent1 = boost::adaptors::transform(g.transitions(),
        std::mem_fn(&RelationshipGraph::transition_t::parent1));
    CHECK_EQUAL_RANGES(test_transition_parent1, expected_transition_parent1);

    auto test_transition_parent2 = boost::adaptors::transform(g.transitions(),
        std::mem_fn(&RelationshipGraph::transition_t::parent2));
    CHECK_EQUAL_RANGES(test_transition_parent2, expected_transition_parent2);

    auto test_transition_length1 = boost::adaptors::transform(g.transitions(),
        std::mem_fn(&RelationshipGraph::transition_t::length1));
    CHECK_EQUAL_RANGES(test_transition_length1, expected_transition_length1);

    auto test_transition_length2 = boost::adaptors::transform(g.transitions(),
        std::mem_fn(&RelationshipGraph::transition_t::length2));
    CHECK_EQUAL_RANGES(test_transition_length2, expected_transition_length2);
}