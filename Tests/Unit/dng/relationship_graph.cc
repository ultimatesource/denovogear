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
#include <utility>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/remove_copy.hpp>
#include <boost/range/irange.hpp>
#include <boost/mem_fn.hpp>

#include "../testing.h"

namespace dng {
struct unittest_dng_relationship_graph {
    GETTERS_FOR_MEMBER_VARIABLES(RelationshipGraph,
        (first_founder)(first_nonfounder)
        (first_somatic)(first_library)(roots)(num_nodes)
        (ploidies)
        (peeling_ops)(peeling_functions_ops)
        (peeling_functions)(peeling_reverse_functions)
        (family_members)
    )
};

template<typename CharType, typename CharTrait>
inline
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const RelationshipGraph::TransitionType& m) {
    o << (int)m;
    return o;
}

namespace peel {
template<typename CharType, typename CharTrait>
inline
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const Op& m) {
    o << (int)m;
    return o;
}
}

}

namespace std {
template<typename CharType, typename CharTrait, typename A, typename B>
inline
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const std::pair<A,B>& m) {
    o << "{" << m.first << ", " << m.second << "}";
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

struct graph_node_t {
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

using expected_graph_nodes_t = std::vector<graph_node_t>;

BOOST_AUTO_TEST_CASE(test_RelationshipGraph_Construct_nodes) {
    const double prec = 2*DBL_EPSILON;

    auto test = [prec](const RelationshipGraph &test, const expected_graph_nodes_t& expected) -> void {
        using boost::adaptors::filtered;
        using boost::adaptors::transformed;

        auto indexes = boost::irange(expected_graph_nodes_t::size_type(0), expected.size());

        // Check Node Structure
        BOOST_CHECK_EQUAL(test.num_nodes(), expected.size());

        auto expected_first_founder = utility::find_position_if(expected, [](const graph_node_t & n) -> bool {
            return (n.vertex_type == GERMLINE && n.transition_type == FOUNDER);
        });
        BOOST_CHECK_EQUAL(u::first_founder(test), expected_first_founder);

        auto expected_first_nonfounder = utility::find_position_if(expected, [](const graph_node_t & n) -> bool {
            return (n.transition_type != FOUNDER);
        });
        BOOST_CHECK_EQUAL(u::first_nonfounder(test), expected_first_nonfounder);

        auto expected_first_somatic = utility::find_position_if(expected, [](const graph_node_t & n) -> bool {
            return (n.vertex_type >= SOMATIC);
        });
        BOOST_CHECK_EQUAL(u::first_somatic(test), expected_first_somatic);

        auto expected_first_library = utility::find_position_if(expected, [](const graph_node_t & n) -> bool {
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
        auto expected_labels = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::label)));
        CHECK_EQUAL_RANGES(test.labels(), expected_labels);

        auto expected_ploidies = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::ploidy)));
        CHECK_EQUAL_RANGES(test.ploidies(), expected_ploidies);

        // Check Library names
        std::vector<std::string> test_library_names = test.library_names();
        boost::sort(test_library_names);
        std::vector<std::string> expected_library_names;
        boost::remove_copy(expected | transformed(boost::mem_fn(&graph_node_t::library_name)),
            std::back_inserter(expected_library_names), "");
        boost::sort(expected_library_names);
        CHECK_EQUAL_RANGES(test_library_names, expected_library_names);

        // Check Transitions
        auto test_transition_types = make_test_range(test.transitions() | transformed(boost::mem_fn(&transition_t::type)));
        auto expected_transition_types = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::transition_type)));
        CHECK_EQUAL_RANGES(test_transition_types, expected_transition_types);

        auto test_parent1 = make_test_range(test.transitions() | transformed(boost::mem_fn(&transition_t::parent1)));
        auto expected_parent1 = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::parent1)));
        CHECK_EQUAL_RANGES(test_parent1, expected_parent1);

        auto test_parent2 = make_test_range(test.transitions() | transformed(boost::mem_fn(&transition_t::parent2)));
        auto expected_parent2 = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::parent2)));
        CHECK_EQUAL_RANGES(test_parent2, expected_parent2);

        auto test_length1 = make_test_range(test.transitions() | transformed(boost::mem_fn(&transition_t::length1)));
        auto expected_length1 = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::length1)));
        CHECK_EQUAL_RANGES(test_length1, expected_length1);

        auto test_length2 = make_test_range(test.transitions() | transformed(boost::mem_fn(&transition_t::length2)));
        auto expected_length2 = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::length2)));
        CHECK_EQUAL_RANGES(test_length2, expected_length2);
    };

    libraries_t quad_libs = {
        {"DadLb", "MomLb","EveLb", "BobLb"},
        {"DadSm", "MomSm", "EveSm", "BobSm"}
    };

    Pedigree quad_ped;
    quad_ped.AddMember("Dad","0","0",Sex::Male,"DadSm");
    quad_ped.AddMember("Mom","0","0",Sex::Female,"MomSm");
    quad_ped.AddMember("Eve","Dad","Mom",Sex::Female,"EveSm");
    quad_ped.AddMember("Bob","Dad","Mom",Sex::Male,"BobSm");

    BOOST_TEST_CONTEXT("graph=quad_graph_autosomal") {
        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::Autosomal, 1e-8, 3e-8, 5e-8, true));

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_graph_nodes_t expected = {
            {{"GL/Dad"}, DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"GL/Mom"}, DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/DadSm/DadLb"}, DIPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"DadLb"}},
            {{"LB/MomSm/MomLb"}, DIPLOID, LIBRARY, !ROOT, PAIR, 1, ls, -1, 0.0f, {"MomLb"}},
            {{"LB/EveSm/EveLb"}, DIPLOID, LIBRARY, !ROOT, TRIO, 0, lsg, 1, lsg, {"EveLb"}},
            {{"LB/BobSm/BobLb"}, DIPLOID, LIBRARY, !ROOT, TRIO, 0, lsg, 1, lsg, {"BobLb"}}
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_xlinked") {
        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::XLinked, 1e-8, 3e-8, 5e-8, true));

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_graph_nodes_t expected = {
            {{"GL/Dad"}, HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"GL/Mom"}, DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/DadSm/DadLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"DadLb"}},
            {{"LB/MomSm/MomLb"}, DIPLOID, LIBRARY, !ROOT, PAIR, 1, ls, -1, 0.0f, {"MomLb"}},
            {{"LB/EveSm/EveLb"}, DIPLOID, LIBRARY, !ROOT, TRIO, 0, lsg, 1, lsg, {"EveLb"}},
            {{"LB/BobSm/BobLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 1, lsg, -1, 0.0f, {"BobLb"}}
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_ylinked") {
        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::YLinked, 1e-8, 3e-8, 5e-8, true));

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_graph_nodes_t expected = {
            {{"GL/Dad"}, HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/DadSm/DadLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"DadLb"}},
            {{"LB/BobSm/BobLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, lsg, -1, 0.0f, {"BobLb"}}
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_zlinked") {
        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::ZLinked, 1e-8, 3e-8, 5e-8, true));

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_graph_nodes_t expected = {
            {{"GL/Dad"}, DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"GL/Mom"}, HAPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/DadSm/DadLb"}, DIPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"DadLb"}},
            {{"LB/MomSm/MomLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 1, ls, -1, 0.0f, {"MomLb"}},
            {{"LB/EveSm/EveLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, lsg, -1, 0.0f, {"EveLb"}},
            {{"LB/BobSm/BobLb"}, DIPLOID, LIBRARY, !ROOT, TRIO, 0, lsg, 1, lsg, {"BobLb"}}
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_wlinked") {
        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::WLinked, 1e-8, 3e-8, 5e-8, true));

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_graph_nodes_t expected = {
            {{"GL/Mom"}, HAPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/MomSm/MomLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"MomLb"}},
            {{"LB/EveSm/EveLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, lsg, -1, 0.0f, {"EveLb"}},
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_maternal") {
        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::Maternal, 1e-8, 3e-8, 5e-8, true));

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_graph_nodes_t expected = {
            {{"GL/Dad"}, HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"GL/Mom"}, HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/DadSm/DadLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"DadLb"}},
            {{"LB/MomSm/MomLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 1, ls, -1, 0.0f, {"MomLb"}},
            {{"LB/EveSm/EveLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 1, lsg, -1, 0.0f, {"EveLb"}},
            {{"LB/BobSm/BobLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 1, lsg, -1, 0.0f, {"BobLb"}}
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_paternal") {
        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::Paternal, 1e-8, 3e-8, 5e-8, true));

        float lsg = 9e-8f; // library-somatic-germline
        float ls  = 8e-8f; // library-somatic

        const expected_graph_nodes_t expected = {
            {{"GL/Dad"}, HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"GL/Mom"}, HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, {""}},
            {{"LB/DadSm/DadLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, ls, -1, 0.0f, {"DadLb"}},
            {{"LB/MomSm/MomLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 1, ls, -1, 0.0f, {"MomLb"}},
            {{"LB/EveSm/EveLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, lsg, -1, 0.0f, {"EveLb"}},
            {{"LB/BobSm/BobLb"}, HAPLOID, LIBRARY, !ROOT, PAIR, 0, lsg, -1, 0.0f, {"BobLb"}}
          };

        test(graph, expected);
    }


    //////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_TEST_CONTEXT("graph=somatic_graph_normalized") {
        libraries_t somatic_libs = {
            {"M1A", "M1B", "M2", "M3A", "M3B"},
            {"M1", "M1", "M2", "M3", "M3"}
        };

        Pedigree somatic_ped;
        somatic_ped.AddMember("M","0","0",Sex::Female,"((M1,M2)M12,M3);");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(somatic_ped, somatic_libs,
            InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/M", DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"SM/unnamed_node_1", DIPLOID, SOMATIC, !ROOT, PAIR, 0, s/3, -1, 0.0f, ""},
            {"SM/M3", DIPLOID, SOMATIC, !ROOT, PAIR, 1, s/3, -1, 0.0f, ""},
            {"SM/M12", DIPLOID, SOMATIC, !ROOT, PAIR, 1, s/3, -1, 0.0f, ""},
            {"SM/M1", DIPLOID, SOMATIC, !ROOT, PAIR, 3, s/3, -1, 0.0f, ""},
            {"LB/M1/M1A", DIPLOID, LIBRARY, !ROOT, PAIR, 4, l, -1, 0.0f, "M1A"},
            {"LB/M1/M1B", DIPLOID, LIBRARY, !ROOT, PAIR, 4, l, -1, 0.0f, "M1B"},            
            {"LB/M2", DIPLOID, LIBRARY, !ROOT, PAIR, 3, l+s/3, -1, 0.0f, "M2"},
            {"LB/M3/M3A", DIPLOID, LIBRARY, !ROOT, PAIR, 2, l, -1, 0.0f, "M3A"},
            {"LB/M3/M3B", DIPLOID, LIBRARY, !ROOT, PAIR, 2, l, -1, 0.0f, "M3B"},
          };

        test(graph, expected);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_TEST_CONTEXT("graph=somatic_graph_unnormalized") {
        libraries_t somatic_libs = {
            {"M1A", "M1B", "M2", "M3A", "M3B"},
            {"M1", "M1", "M2", "M3", "M3"}
        };

        Pedigree somatic_ped;
        somatic_ped.AddMember("M","0","0",Sex::Female,"((M1,M2)M12,M3);");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(somatic_ped, somatic_libs,
            InheritanceModel::Autosomal, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/M", DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"SM/unnamed_node_1", DIPLOID, SOMATIC, !ROOT, PAIR, 0, s, -1, 0.0f, ""},
            {"SM/M3", DIPLOID, SOMATIC, !ROOT, PAIR, 1, s, -1, 0.0f, ""},
            {"SM/M12", DIPLOID, SOMATIC, !ROOT, PAIR, 1, s, -1, 0.0f, ""},
            {"SM/M1", DIPLOID, SOMATIC, !ROOT, PAIR, 3, s, -1, 0.0f, ""},
            {"LB/M1/M1A", DIPLOID, LIBRARY, !ROOT, PAIR, 4, l, -1, 0.0f, "M1A"},
            {"LB/M1/M1B", DIPLOID, LIBRARY, !ROOT, PAIR, 4, l, -1, 0.0f, "M1B"},            
            {"LB/M2", DIPLOID, LIBRARY, !ROOT, PAIR, 3, l+s, -1, 0.0f, "M2"},
            {"LB/M3/M3A", DIPLOID, LIBRARY, !ROOT, PAIR, 2, l, -1, 0.0f, "M3A"},
            {"LB/M3/M3B", DIPLOID, LIBRARY, !ROOT, PAIR, 2, l, -1, 0.0f, "M3B"},
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=trio_graph") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember("Mom","0","0",Sex::Female,"");
        ped.AddMember("Dad","0","0",Sex::Male,"");
        ped.AddMember("Eve","Dad","Mom",Sex::Female,"");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Mom", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Dad", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, TRIO, 1, s+g+l, 0, s+g+l, "Eve"},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Mom"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=trio_graph_no_mom") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember("Dad","0","0",Sex::Male,"");
        ped.AddMember("Eve","Dad","Mom",Sex::Female,"");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/unknown_mom_of_Eve", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, TRIO, 0, s+g+l, 1, s+g+l, "Eve"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=trio_graph_no_mom_lib") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember("Eve","Dad","Mom2",Sex::Female,"");
        ped.AddMember("Dad","0","0",Sex::Male,"");
        ped.AddMember("Mom2","0","0",Sex::Female,"");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom2", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, TRIO, 0, s+g+l, 1, s+g+l, "Eve"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);        
    }


    BOOST_TEST_CONTEXT("graph=trio_graph_no_dad") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember("Mom","0","0",Sex::Female,"");
        ped.AddMember("Eve","Dad","Mom",Sex::Female,"");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Mom", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/unknown_dad_of_Eve", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, TRIO, 1, s+g+l, 0, s+g+l, "Eve"},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Mom"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=trio_graph_no_dad_lib") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember("Mom","0","0",Sex::Female,"");
        ped.AddMember("Eve","Dad2","Mom",Sex::Female,"");
        ped.AddMember("Dad2","0","0",Sex::Male,"");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Mom", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Dad2", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, TRIO, 1, s+g+l, 0, s+g+l, "Eve"},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Mom"},
          };

        test(graph, expected);        
    }


    BOOST_TEST_CONTEXT("graph=trio_graph_no_child") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember("Mom","0","0",Sex::Female,"");
        ped.AddMember("Dad","0","0",Sex::Male,"");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Mom", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Mom"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=trio_graph_no_child_lib") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember("Mom","0","0",Sex::Female,"");
        ped.AddMember("Dad","0","0",Sex::Male,"");
        ped.AddMember("Eve2","Dad","Mom",Sex::Female,"");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Mom", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Mom"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=halfsibs_graph") {
        libraries_t libs = {
            {"Dad", "Mom1", "Mom2", "Eve1", "Bob2"},
            {"Dad", "Mom1", "Mom2", "Eve1", "Bob2"}
        };

        Pedigree ped;
        ped.AddMember("Dad","0","0",Sex::Male,"");
        ped.AddMember("Mom1","0","0",Sex::Female,"");
        ped.AddMember("Mom2","0","0",Sex::Female,"");
        ped.AddMember("Eve1","Dad","Mom1",Sex::Female,"");
        ped.AddMember("Bob2","Dad","Mom2",Sex::Male,"");

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom1", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom2", DIPLOID, GERMLINE,  !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
            {"LB/Mom1", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom1"},
            {"LB/Mom2", DIPLOID, LIBRARY,  !ROOT, PAIR, 2, s+l, -1, 0.0f, "Mom2"},
            {"LB/Eve1", DIPLOID, LIBRARY,  !ROOT, TRIO, 0, g+s+l, 1, g+s+l, "Eve1"},
            {"LB/Bob2", DIPLOID, LIBRARY,  !ROOT, TRIO, 0, g+s+l, 2, g+s+l, "Bob2"},
          };

        test(graph, expected);
    }
}

BOOST_AUTO_TEST_CASE(test_RelationshipGraph_CreateWorkspace) {
    // Test CreateWorkspace 

    using std::make_pair;
    RelationshipGraph graph;
    u::first_founder(graph) = 0;
    u::first_nonfounder(graph) = 2;
    u::first_somatic(graph) = 4; 
    u::first_library(graph) = 8;
    u::num_nodes(graph) = 12;
    u::ploidies(graph).assign(12,2);

    auto workspace = graph.CreateWorkspace();

    using p = decltype(workspace.founder_nodes);

    BOOST_CHECK_EQUAL(workspace.num_nodes, 12);
    BOOST_CHECK_EQUAL(workspace.founder_nodes,  (p{0,2}));
    BOOST_CHECK_EQUAL(workspace.germline_nodes, (p{0,4}));
    BOOST_CHECK_EQUAL(workspace.somatic_nodes,  (p{4,8}));
    BOOST_CHECK_EQUAL(workspace.library_nodes,  (p{8,12}));
    
    std::vector<int> expected_ploidies(12,2);
    CHECK_EQUAL_RANGES(workspace.ploidies, expected_ploidies);

    BOOST_CHECK_EQUAL(workspace.upper.size(), 12);
    BOOST_CHECK_EQUAL(workspace.lower.size(), 12);
    BOOST_CHECK_EQUAL(workspace.super.size(), 12);
    BOOST_CHECK_EQUAL(workspace.dirty_lower, false);
}

struct peeling_node_t {
    peel::Op op;
    peel::family_members_t family;
};

using expected_peeling_nodes_t = std::vector<peeling_node_t>;

BOOST_AUTO_TEST_CASE(test_RelationshipGraph_Construct_peeler) {
    const double prec = 2*DBL_EPSILON;

    auto test = [prec](const RelationshipGraph &test, const expected_peeling_nodes_t& expected) -> void {
        using boost::adaptors::filtered;
        using boost::adaptors::transformed;
        using Op = peel::Op;
        using function_t = peel::function_t;

        auto test_fastops = u::peeling_functions_ops(test);
        auto expected_fastops = make_test_range(expected | transformed(boost::mem_fn(&peeling_node_t::op)));
        CHECK_EQUAL_RANGES(test_fastops, expected_fastops);

        auto test_origops = u::peeling_ops(test);
        auto expected_origops = make_test_range(expected | transformed([](const peeling_node_t &n) -> Op {
             return ((int)n.op < (int)Op::UPFAST) ? n.op : (Op)((int)n.op - (int)Op::UPFAST); 
        }));
        CHECK_EQUAL_RANGES(test_origops, expected_origops);

        auto test_forward = u::peeling_functions(test);
        auto expected_forward = make_test_range(expected | transformed([](const peeling_node_t &n) -> function_t {
            switch(n.op) {
            case Op::UP:
                return &peel::up;
            case Op::DOWN:
                return &peel::down;
            case Op::TOFATHER:
                return &peel::to_father;
            case Op::TOMOTHER:
                return &peel::to_mother;
            case Op::TOCHILD:
                return &peel::to_child;
            case Op::UPFAST:
                return &peel::up_fast;
            case Op::DOWNFAST:
                return &peel::down_fast;
            case Op::TOFATHERFAST:
                return &peel::to_father_fast;
            case Op::TOMOTHERFAST:
                return &peel::to_mother_fast;
            case Op::TOCHILDFAST:
                return &peel::to_child_fast;
            default:
                return NULL;
            };
        }));
        CHECK_EQUAL_RANGES(test_forward, expected_forward);

        auto test_reverse = u::peeling_reverse_functions(test);
        auto expected_reverse = make_test_range(expected | transformed([](const peeling_node_t &n) -> function_t {
            switch(n.op) {
            case Op::UP:
                return &peel::up_reverse;
            case Op::DOWN:
                return &peel::down_reverse;
            case Op::TOFATHER:
                return &peel::to_father_reverse;
            case Op::TOMOTHER:
                return &peel::to_mother_reverse;
            case Op::TOCHILD:
                return &peel::to_child_reverse;
            case Op::UPFAST:
                return &peel::up_reverse;
            case Op::DOWNFAST:
                return &peel::down_reverse;
            case Op::TOFATHERFAST:
                return &peel::to_father_reverse;
            case Op::TOMOTHERFAST:
                return &peel::to_mother_reverse;
            case Op::TOCHILDFAST:
                return &peel::to_child_reverse;
            default:
                return NULL;
            };
        }));
        CHECK_EQUAL_RANGES(test_reverse, expected_reverse);

        auto test_families = u::family_members(test); 
        auto expected_families = expected | transformed(boost::mem_fn(&peeling_node_t::family));
        BOOST_REQUIRE_EQUAL(test_families.size(), expected_families.size());
        for(int family_number=0; family_number < test_families.size(); ++family_number){
            BOOST_TEST_INFO("family_number=" << family_number);
            auto test_family = test_families[family_number];
            auto expected_family = expected_families[family_number];
            CHECK_EQUAL_RANGES(test_family, expected_family);
        }
    };

    using Op = peel::Op;

    constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

    libraries_t quad_libs = {
        {"DadLb", "MomLb","EveLb", "BobLb"},
        {"DadSm", "MomSm", "EveSm", "BobSm"}
    };

    Pedigree quad_ped;
    quad_ped.AddMember("Dad","0","0",Sex::Male,"DadSm");
    quad_ped.AddMember("Mom","0","0",Sex::Female,"MomSm");
    quad_ped.AddMember("Eve","Dad","Mom",Sex::Female,"EveSm");
    quad_ped.AddMember("Bob","Dad","Mom",Sex::Male,"BobSm");

    expected_peeling_nodes_t expected = {
        {Op::UPFAST, {1,3}},
        {Op::TOFATHERFAST, {0,1,5,4}},
        {Op::UP, {0,2}},
    };

    RelationshipGraph graph;
    BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
        InheritanceModel::Autosomal, g, s, l, true));
    test(graph, expected);
}
