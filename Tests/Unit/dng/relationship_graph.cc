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
#include <dng/utility.h>

#include <vector>
#include <functional>
#include <iterator>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/remove_copy.hpp>
#include <boost/range/irange.hpp>

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

struct node_t {
    std::string label;
    int ploidy;
    int vertex_type;
    bool is_root;
    TransitionType transition_type;
    int parent1; // dad or source
    float length1;
    int parent2; // mom or null
    float length2;
    std::string library_name;
};

constexpr int HAPLOID = 1;
constexpr int DIPLOID = 2;

constexpr int GERMLINE = 1;
constexpr int SOMATIC = 2;
constexpr int LIBRARY = 3;

constexpr auto FOUNDER = TransitionType::Founder;
constexpr auto PAIR = TransitionType::Pair;
constexpr auto TRIO = TransitionType::Trio;

constexpr bool ROOT = true;

using expected_nodes_t = std::vector<node_t>;

BOOST_AUTO_TEST_CASE(test_RelationshipGraph_Construct) {
    const double prec = 2*DBL_EPSILON;

    auto test = [prec](const RelationshipGraph &test, const expected_nodes_t& expected) -> void {
        using boost::adaptors::filtered;
        using boost::adaptors::transformed;

        auto indexes = boost::irange(expected_nodes_t::size_type(0), expected.size());

        // Check Node Structure
        BOOST_CHECK_EQUAL(test.num_nodes(), expected.size());

        auto expected_first_founder = utility::find_position_if(expected, [](const node_t & n) -> bool {
            return (n.vertex_type == GERMLINE && n.transition_type == FOUNDER);
        });
        BOOST_CHECK_EQUAL(u::first_founder(test), expected_first_founder);

        auto expected_first_nonfounder = utility::find_position_if(expected, [](const node_t & n) -> bool {
            return (n.transition_type != FOUNDER);
        });
        BOOST_CHECK_EQUAL(u::first_nonfounder(test), expected_first_nonfounder);

        auto expected_first_somatic = utility::find_position_if(expected, [](const node_t & n) -> bool {
            return (n.vertex_type >= SOMATIC);
        });
        BOOST_CHECK_EQUAL(u::first_somatic(test), expected_first_somatic);

        auto expected_first_library = utility::find_position_if(expected, [](const node_t & n) -> bool {
            return (n.vertex_type >= LIBRARY);
        });
        BOOST_CHECK_EQUAL(u::first_library(test), expected_first_library);

        BOOST_CHECK_EQUAL(test.library_nodes().first, expected_first_library);
        BOOST_CHECK_EQUAL(test.library_nodes().second, expected.size());

        // BOOST_TEST(a == b) doesn't like its second operator to be a boost::range
        // Wrapping the range hides the type and allows the BOOST_TEST operator to be found 
        auto expected_roots_ = indexes | filtered([&](size_t i) -> bool {
            return (expected[i].is_root == true);
        });
        auto expected_roots = make_test_range(expected_roots_); 
        CHECK_EQUAL_RANGES(u::roots(test), expected_roots);

        // Check Labels and Ploidies
        auto expected_labels = make_test_range(expected | transformed(std::mem_fn(&node_t::label)));
        CHECK_EQUAL_RANGES(test.labels(), expected_labels);

        auto expected_ploidies = make_test_range(expected | transformed(std::mem_fn(&node_t::ploidy)));
        CHECK_EQUAL_RANGES(test.ploidies(), expected_ploidies);

        // Check Library names
        std::vector<std::string> test_library_names = test.library_names();
        boost::sort(test_library_names);
        std::vector<std::string> expected_library_names;
        boost::remove_copy(expected | transformed(std::mem_fn(&node_t::library_name)),
            std::back_inserter(expected_library_names), "");
        boost::sort(expected_library_names);
        CHECK_EQUAL_RANGES(test_library_names, expected_library_names);

        // Check Transitions
        auto test_transition_types = make_test_range(test.transitions() | transformed(std::mem_fn(&transition_t::type)));
        auto expected_transition_types = make_test_range(expected | transformed(std::mem_fn(&node_t::transition_type)));
        CHECK_EQUAL_RANGES(test_transition_types, expected_transition_types);

        auto test_parent1 = make_test_range(test.transitions() | transformed(std::mem_fn(&transition_t::parent1)));
        auto expected_parent1 = make_test_range(expected | transformed(std::mem_fn(&node_t::parent1)));
        CHECK_EQUAL_RANGES(test_parent1, expected_parent1);

        auto test_parent2 = make_test_range(test.transitions() | transformed(std::mem_fn(&transition_t::parent2)));
        auto expected_parent2 = make_test_range(expected | transformed(std::mem_fn(&node_t::parent2)));
        CHECK_EQUAL_RANGES(test_parent2, expected_parent2);

        auto test_length1 = make_test_range(test.transitions() | transformed(std::mem_fn(&transition_t::length1)));
        auto expected_length1 = make_test_range(expected | transformed(std::mem_fn(&node_t::length1)));
        CHECK_EQUAL_RANGES(test_length1, expected_length1);

        auto test_length2 = make_test_range(test.transitions() | transformed(std::mem_fn(&transition_t::length2)));
        auto expected_length2 = make_test_range(expected | transformed(std::mem_fn(&node_t::length2)));
        CHECK_EQUAL_RANGES(test_length2, expected_length2);
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
        //Construct Graph
        RelationshipGraph graph;
        graph.Construct(quad_ped, libs, InheritanceModel::Autosomal, 1e-8, 3e-8, 5e-8);

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_nodes_t expected = {
            {{"GL/Dad"}, DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"GL/Mom"}, DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/DadSm/DadLb"}, DIPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"DadLb"}},
            {{"LB/MomSm/MomLb"}, DIPLOID, LIBRARY, !ROOT, PAIR, 1, ls, -1, 0.0f, {"MomLb"}},
            {{"LB/EveSm/EveLb"}, DIPLOID, LIBRARY, !ROOT, TRIO, 0, lsg, 1, lsg, {"EveLb"}},
            {{"LB/BobSm/BobLb"}, DIPLOID, LIBRARY, !ROOT, TRIO, 0, lsg, 1, lsg, {"BobLb"}}
          };

        BOOST_TEST_CONTEXT("graph=quad_graph_autosomal") {
            test(graph, expected);
        }
    }

    {
        //Construct Graph
        RelationshipGraph graph;
        graph.Construct(quad_ped, libs, InheritanceModel::XLinked, 1e-8, 3e-8, 5e-8);

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_nodes_t expected = {
            {{"GL/Dad"}, HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"GL/Mom"}, DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/DadSm/DadLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"DadLb"}},
            {{"LB/MomSm/MomLb"}, DIPLOID, LIBRARY, !ROOT, PAIR, 1, ls, -1, 0.0f, {"MomLb"}},
            {{"LB/EveSm/EveLb"}, DIPLOID, LIBRARY, !ROOT, TRIO, 0, lsg, 1, lsg, {"EveLb"}},
            {{"LB/BobSm/BobLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 1, lsg, -1, 0.0f, {"BobLb"}}
          };

        BOOST_TEST_CONTEXT("graph=quad_graph_xlinked") {
            test(graph, expected);
        }
    }

    {
        //Construct Graph
        RelationshipGraph graph;
        graph.Construct(quad_ped, libs, InheritanceModel::YLinked, 1e-8, 3e-8, 5e-8);

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_nodes_t expected = {
            {{"GL/Dad"}, HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/DadSm/DadLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"DadLb"}},
            {{"LB/BobSm/BobLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, lsg, -1, 0.0f, {"BobLb"}}
          };

        BOOST_TEST_CONTEXT("graph=quad_graph_ylinked") {
            test(graph, expected);
        }
    }

    libraries_t somatic_libs = {
        {"M1A", "M1B", "M2A", "M3A", "M3B"},
        {"M1", "M1", "M2", "M3", "M3"}
    };

    Pedigree somatic_ped;
    somatic_ped.AddMember("M","0","0",Sex::Female,"((M1,M2)M12,M3);");

    {
        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        graph.Construct(somatic_ped, somatic_libs, InheritanceModel::Autosomal, g, s, l);


        const expected_nodes_t expected = {
            {"GL/M", DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"SM/unnamed_node_1", DIPLOID, SOMATIC, !ROOT, PAIR, 0, s, -1, 0.0f, ""},
            {"SM/M3", DIPLOID, SOMATIC, !ROOT, PAIR, 1, s, -1, 0.0f, ""},
            {"SM/M12", DIPLOID, SOMATIC, !ROOT, PAIR, 1, s, -1, 0.0f, ""},
            {"SM/M1", DIPLOID, SOMATIC, !ROOT, PAIR, 3, s, -1, 0.0f, ""},
            {"LB/M1/M1A", DIPLOID, LIBRARY, !ROOT, PAIR, 4, l, -1, 0.0f, "M1A"},
            {"LB/M1/M1B", DIPLOID, LIBRARY, !ROOT, PAIR, 4, l, -1, 0.0f, "M1B"},            
            {"LB/M2/M2A", DIPLOID, LIBRARY, !ROOT, PAIR, 3, l+s, -1, 0.0f, "M2A"},
            {"LB/M3/M3A", DIPLOID, LIBRARY, !ROOT, PAIR, 2, l, -1, 0.0f, "M3A"},
            {"LB/M3/M3B", DIPLOID, LIBRARY, !ROOT, PAIR, 2, l, -1, 0.0f, "M3B"},
          };

        BOOST_TEST_CONTEXT("graph=somatic_graph") {
            test(graph, expected);
        }
    }


}
