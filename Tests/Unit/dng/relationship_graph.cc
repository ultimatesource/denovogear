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
#include <dng/mutation.h>

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

#include <boost/range/algorithm/generate.hpp>
#include <boost/range/algorithm/copy.hpp>

#include "../testing.h"
#include "../xorshift64.h"

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
    o << static_cast<int>(m);
    return o;
}

namespace peel {
template<typename CharType, typename CharTrait>
inline
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const Op& m) {
    o << static_cast<int>(m);
    return o;
}

template<typename CharType, typename CharTrait>
inline
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, function_t m) {
    o << reinterpret_cast<void*>(m);
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

const int NUM_TEST = 25;
int g_seed_counter = 0;

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

        std::vector<std::string> test_labels_pos;
        for(int i=0; i< test.num_nodes(); ++i) {
            test_labels_pos.push_back(test.label(i));
        }
        CHECK_EQUAL_RANGES(test_labels_pos, expected_labels);

        auto expected_ploidies = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::ploidy)));
        CHECK_EQUAL_RANGES(test.ploidies(), expected_ploidies);

        std::vector<int> test_ploidies_pos;
        for(int i=0; i< test.num_nodes(); ++i) {
            test_ploidies_pos.push_back(test.ploidy(i));
        }
        CHECK_EQUAL_RANGES(test_ploidies_pos, expected_ploidies);

        // Check Library names
        std::vector<std::string> test_library_names = test.library_names();
        boost::sort(test_library_names);
        std::vector<std::string> expected_library_names;
        boost::remove_copy(expected | transformed(boost::mem_fn(&graph_node_t::library_name)),
            std::back_inserter(expected_library_names), "");
        boost::sort(expected_library_names);
        CHECK_EQUAL_RANGES(test_library_names, expected_library_names);

        // Check Transitions
        auto transitions_test = [&](const std::vector<transition_t> & test_transitions,
                const char * msg) {
            BOOST_TEST_CONTEXT(msg) {
            
            auto test_transition_types = make_test_range(test_transitions | transformed(boost::mem_fn(&transition_t::type)));
            auto expected_transition_types = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::transition_type)));
            CHECK_EQUAL_RANGES(test_transition_types, expected_transition_types);

            auto test_parent1 = make_test_range(test_transitions | transformed(boost::mem_fn(&transition_t::parent1)));
            auto expected_parent1 = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::parent1)));
            CHECK_EQUAL_RANGES(test_parent1, expected_parent1);

            auto test_parent2 = make_test_range(test_transitions | transformed(boost::mem_fn(&transition_t::parent2)));
            auto expected_parent2 = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::parent2)));
            CHECK_EQUAL_RANGES(test_parent2, expected_parent2);

            auto test_length1 = make_test_range(test_transitions | transformed(boost::mem_fn(&transition_t::length1)));
            auto expected_length1 = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::length1)));
            CHECK_EQUAL_RANGES(test_length1, expected_length1);

            auto test_length2 = make_test_range(test_transitions | transformed(boost::mem_fn(&transition_t::length2)));
            auto expected_length2 = make_test_range(expected | transformed(boost::mem_fn(&graph_node_t::length2)));
            CHECK_EQUAL_RANGES(test_length2, expected_length2);
        }};
        std::vector<transition_t> test_transitions_pos;
        for(int i=0; i< test.num_nodes(); ++i) {
            test_transitions_pos.push_back(test.transition(i));
        }

        transitions_test(test.transitions(), "test.transitions()");
        transitions_test(test_transitions_pos, "test_transitions_pos");
    };

    libraries_t quad_libs = {
        {"DadLb", "MomLb","EveLb", "BobLb"},
        {"DadSm", "MomSm", "EveSm", "BobSm"}
    };

    Pedigree quad_ped;
    quad_ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"DadSm"}});
    quad_ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"MomSm"}});
    quad_ped.AddMember({"Eve",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Female,{"EveSm"}});
    quad_ped.AddMember({"Bob",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Male,{"BobSm"}});

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
            {{"LB/BobSm/BobLb"}, DIPLOID, LIBRARY, !ROOT, TRIO, 0, lsg, 1, lsg, {"BobLb"}},
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
        somatic_ped.AddMember({"M",{},{},{},{},{},Sex::Autosomal,{"((M1,M2)M12,M3);"}});

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
        somatic_ped.AddMember({"M",{},{},{},{},{},Sex::Autosomal,{"((M1,M2)M12,M3);"}});

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

    BOOST_TEST_CONTEXT("graph=somatic_graph_star") {
        libraries_t somatic_libs = {
            {"M1A", "M1B", "M2", "M3A", "M3B"},
            {"M1", "M1", "M2", "M3", "M3"}
        };

        Pedigree somatic_ped;
        somatic_ped.AddMember({"M",{},{},{},{},{},Sex::Autosomal,{"M1","M2","M3:0.0"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(somatic_ped, somatic_libs,
            InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/M",      DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"SM/M1",     DIPLOID, SOMATIC, !ROOT, PAIR, 0, s, -1, 0.0f, ""},
            {"SM/M3",     DIPLOID, SOMATIC, !ROOT, PAIR, 0, 0.0, -1, 0.0f, ""},
            {"LB/M1/M1A", DIPLOID, LIBRARY, !ROOT, PAIR, 1, l, -1, 0.0f, "M1A"},
            {"LB/M1/M1B", DIPLOID, LIBRARY, !ROOT, PAIR, 1, l, -1, 0.0f, "M1B"},            
            {"LB/M2",     DIPLOID, LIBRARY, !ROOT, PAIR, 0, l+s, -1, 0.0f, "M2"},
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
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"Eve",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Female,{"Eve"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, TRIO, 0, s+g+l, 1, s+g+l, "Eve"},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
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
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"MomSm"}});
        ped.AddMember({"Eve",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Female,{"Eve"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, TRIO, 0, s+g+l, 1, s+g+l, "Eve"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=trio_graph_no_dad_lib") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"DadSm"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"Eve",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Female,{"Eve"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, TRIO, 0, s+g+l, 1, s+g+l, "Eve"},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
          };

        test(graph, expected);
    }


    BOOST_TEST_CONTEXT("graph=trio_graph_no_child") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=trio_graph_no_child_lib") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"Eve",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Female,{"EveSm"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=trio_graph_child_is_founder") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"Eve",{"founder"},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Female,{"Eve"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Eve", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, PAIR, 2, s+l, -1, 0.0f, "Eve"},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);        
    }


    BOOST_TEST_CONTEXT("graph=halfsibs_graph") {
        libraries_t libs = {
            {"Dad", "Mom1", "Mom2", "Eve1", "Bob2"},
            {"Dad", "Mom1", "Mom2", "Eve1", "Bob2"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom1",{},{},{},{},{},Sex::Female,{"Mom1"}});
        ped.AddMember({"Mom2",{},{},{},{},{},Sex::Female,{"Mom2"}});
        ped.AddMember({"Eve1",{},std::string{"Dad"},{},std::string{"Mom1"},{},Sex::Female,{"Eve1"}});
        ped.AddMember({"Bob2",{},std::string{"Dad"},{},std::string{"Mom2"},{},Sex::Male,{"Bob2"}});

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

    BOOST_TEST_CONTEXT("graph=m12") {
/*
1-2    3-4
 |      |
 7------8    5-6
   |  |       |
   9  10-----11
          |
          12
*/
        libraries_t libs = {
            {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"},
            {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"}
        };

        Pedigree ped;
        ped.AddMember({"1",{},{},{},{},{},Sex::Male,{"1"}});
        ped.AddMember({"2",{},{},{},{},{},Sex::Female,{"2"}});
        ped.AddMember({"3",{},{},{},{},{},Sex::Male,{"3"}});
        ped.AddMember({"4",{},{},{},{},{},Sex::Female,{"4"}});
        ped.AddMember({"5",{},{},{},{},{},Sex::Male,{"5"}});
        ped.AddMember({"6",{},{},{},{},{},Sex::Female,{"6"}});
        ped.AddMember({"7",{},std::string{"1"},{},std::string{"2"},{},Sex::Male,{"7"}});
        ped.AddMember({"8",{},std::string{"3"},{},std::string{"4"},{},Sex::Female,{"8"}});
        ped.AddMember({"9",{},std::string{"7"},{},std::string{"8"},{},Sex::Female,{"9"}});
        ped.AddMember({"10",{},std::string{"7"},{},std::string{"8"},{},Sex::Female,{"10"}});
        ped.AddMember({"11",{},std::string{"5"},{},std::string{"6"},{},Sex::Male,{"11"}});
        ped.AddMember({"12",{},std::string{"11"},{},std::string{"10"},{},Sex::Male,{"12"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/1", DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/2", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/3", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/4", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/5", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/6", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            
            {"GL/7", DIPLOID, GERMLINE, !ROOT, TRIO, 0, g, 1, g, ""},
            {"GL/8", DIPLOID, GERMLINE, !ROOT, TRIO, 2, g, 3, g, ""},
            {"GL/11", DIPLOID, GERMLINE, !ROOT, TRIO, 4, g, 5, g, ""},
            {"GL/10", DIPLOID, GERMLINE, !ROOT, TRIO, 6, g, 7, g, ""},

            {"LB/1", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "1"},
            {"LB/2", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "2"},
            {"LB/3", DIPLOID, LIBRARY,  !ROOT, PAIR, 2, s+l, -1, 0.0f, "3"},
            {"LB/4", DIPLOID, LIBRARY,  !ROOT, PAIR, 3, s+l, -1, 0.0f, "4"},
            {"LB/5", DIPLOID, LIBRARY,  !ROOT, PAIR, 4, s+l, -1, 0.0f, "5"},
            {"LB/6", DIPLOID, LIBRARY,  !ROOT, PAIR, 5, s+l, -1, 0.0f, "6"},
            {"LB/7", DIPLOID, LIBRARY,  !ROOT, PAIR, 6, s+l, -1, 0.0f, "7"},
            {"LB/8", DIPLOID, LIBRARY,  !ROOT, PAIR, 7, s+l, -1, 0.0f, "8"},
            {"LB/9", DIPLOID, LIBRARY,  !ROOT, TRIO, 6, g+s+l, 7, g+s+l, "9"},
            {"LB/10", DIPLOID, LIBRARY,  !ROOT, PAIR, 9, s+l, -1, 0.0f, "10"},
            {"LB/11", DIPLOID, LIBRARY,  !ROOT, PAIR, 8, s+l, -1, 0.0f, "11"},
            {"LB/12", DIPLOID, LIBRARY,  !ROOT, TRIO, 8, g+s+l, 9, g+s+l, "12"},
          };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=trio_graph_branch_lengths") {
        libraries_t libs = {
            {"Eve", "Mom", "Dad"},
            {"Eve", "Mom", "Dad"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad:0.1"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"Eve",{},std::string{"Dad"},2,std::string{"Mom"},3,Sex::Female,{"Eve"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/Dad", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom", DIPLOID, GERMLINE, !ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Eve", DIPLOID, LIBRARY,  !ROOT, TRIO, 0, s+2*g+l, 1, s+3*g+l, "Eve"},
            {"LB/Mom", DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/Dad", DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s*0.1+l, -1, 0.0f, "Dad"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=gamete_graph_autosomal") {
        libraries_t libs = {
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"},
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"SpermEve",{"gamete"},std::string{"Dad"},{},{},{},Sex::Female,{"SpermEve"}});
        ped.AddMember({"EggEve",  {"gamete"},{},{},std::string{"Mom"},{},Sex::Female,{"EggEve"}});
        ped.AddMember({"SpermBob",{"haploid"},std::string{"Dad"},{},{},{},Sex::Male,{"SpermBob"}});
        ped.AddMember({"EggBob",  {"haploid"},{},{},std::string{"Mom"},{},Sex::Male,{"EggBob"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/Dad",      DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom",      DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Dad",      DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
            {"LB/Mom",      DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/SpermEve", HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+g+l, -1, 0.0f, "SpermEve"},
            {"LB/EggEve",   HAPLOID, LIBRARY,  !ROOT, PAIR, 1, s+g+l, -1, 0.0f, "EggEve"},
            {"LB/SpermBob", HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+g+l, -1, 0.0f, "SpermBob"},
            {"LB/EggBob",   HAPLOID, LIBRARY,  !ROOT, PAIR, 1, s+g+l, -1, 0.0f, "EggBob"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=gamete_graph_ylinked") {
        libraries_t libs = {
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"},
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"SpermEve",{"gamete"},std::string{"Dad"},{},{},{},Sex::Female,{"SpermEve"}});
        ped.AddMember({"EggEve",  {"gamete"},{},{},std::string{"Mom"},{},Sex::Female,{"EggEve"}});
        ped.AddMember({"SpermBob",{"haploid"},std::string{"Dad"},{},{},{},Sex::Male,{"SpermBob"}});
        // Under xy-linkage, the following is invalid and will be removed because 'male'-gametes must come from dads.
        ped.AddMember({"EggBob",  {"haploid"},{},{},std::string{"Mom"},{},Sex::Male,{"EggBob"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::YLinked, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/Dad",      HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Dad",      HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
            {"LB/SpermBob", HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+g+l, -1, 0.0f, "SpermBob"},
          };

        test(graph, expected);        
    }


    BOOST_TEST_CONTEXT("graph=gamete_graph_xlinked") {
        libraries_t libs = {
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"},
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"SpermEve",{"gamete"},std::string{"Dad"},{},{},{},Sex::Female,{"SpermEve"}});
        ped.AddMember({"EggEve",  {"gamete"},{},{},std::string{"Mom"},{},Sex::Female,{"EggEve"}});
        ped.AddMember({"SpermBob",{"haploid"},std::string{"Dad"},{},{},{},Sex::Male,{"SpermBob"}});
        // Under xy-linkage, the following is invalid and will be removed because 'male'-gametes must come from dads.
        ped.AddMember({"EggBob",  {"haploid"},{},{},std::string{"Mom"},{},Sex::Male,{"EggBob"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::XLinked, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/Dad",      HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom",      DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Dad",      HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
            {"LB/Mom",      DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/SpermEve", HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+g+l, -1, 0.0f, "SpermEve"},
            {"LB/EggEve",   HAPLOID, LIBRARY,  !ROOT, PAIR, 1, s+g+l, -1, 0.0f, "EggEve"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=gamete_graph_wlinked") {
        libraries_t libs = {
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"},
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"SpermEve",{"gamete"},std::string{"Dad"},{},{},{},Sex::Female,{"SpermEve"}});
        ped.AddMember({"EggEve",  {"gamete"},{},{},std::string{"Mom"},{},Sex::Female,{"EggEve"}});
        ped.AddMember({"SpermBob",{"haploid"},std::string{"Dad"},{},{},{},Sex::Male,{"SpermBob"}});
        ped.AddMember({"EggBob",  {"haploid"},{},{},std::string{"Mom"},{},Sex::Male,{"EggBob"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::WLinked, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/Mom",      HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Mom",      HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Mom"},
            {"LB/EggEve",   HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+g+l, -1, 0.0f, "EggEve"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=gamete_graph_zlinked") {
        libraries_t libs = {
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"},
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"SpermEve",{"gamete"},std::string{"Dad"},{},{},{},Sex::Female,{"SpermEve"}});
        ped.AddMember({"EggEve",  {"gamete"},{},{},std::string{"Mom"},{},Sex::Female,{"EggEve"}});
        ped.AddMember({"SpermBob",{"haploid"},std::string{"Dad"},{},{},{},Sex::Male,{"SpermBob"}});
        ped.AddMember({"EggBob",  {"haploid"},{},{},std::string{"Mom"},{},Sex::Male,{"EggBob"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::ZLinked, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/Dad",      DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom",      HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Dad",      DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
            {"LB/Mom",      HAPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/SpermBob", HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+g+l, -1, 0.0f, "SpermBob"},
            {"LB/EggBob",   HAPLOID, LIBRARY,  !ROOT, PAIR, 1, s+g+l, -1, 0.0f, "EggBob"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=gamete_graph_maternal") {
        libraries_t libs = {
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"},
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"SpermEve",{"gamete"},std::string{"Dad"},{},{},{},Sex::Female,{"SpermEve"}});
        ped.AddMember({"EggEve",  {"gamete"},{},{},std::string{"Mom"},{},Sex::Female,{"EggEve"}});
        ped.AddMember({"SpermBob",{"haploid"},std::string{"Dad"},{},{},{},Sex::Male,{"SpermBob"}});
        ped.AddMember({"EggBob",  {"haploid"},{},{},std::string{"Mom"},{},Sex::Male,{"EggBob"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Maternal, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/Dad",      HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom",      HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Dad",      HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
            {"LB/Mom",      HAPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/EggEve",   HAPLOID, LIBRARY,  !ROOT, PAIR, 1, s+g+l, -1, 0.0f, "EggEve"},
            {"LB/EggBob",   HAPLOID, LIBRARY,  !ROOT, PAIR, 1, s+g+l, -1, 0.0f, "EggBob"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=gamete_graph_paternal") {
        libraries_t libs = {
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"},
            {"Dad", "Mom", "SpermEve", "EggEve", "SpermBob", "EggBob"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"SpermEve",{"gamete"},std::string{"Dad"},{},{},{},Sex::Female,{"SpermEve"}});
        ped.AddMember({"EggEve",  {"gamete"},{},{},std::string{"Mom"},{},Sex::Female,{"EggEve"}});
        ped.AddMember({"SpermBob",{"haploid"},std::string{"Dad"},{},{},{},Sex::Male,{"SpermBob"}});
        ped.AddMember({"EggBob",  {"haploid"},{},{},std::string{"Mom"},{},Sex::Male,{"EggBob"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Paternal, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/Dad",      HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom",      HAPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Dad",      HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
            {"LB/Mom",      HAPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/SpermEve", HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+g+l, -1, 0.0f, "SpermEve"},
            {"LB/SpermBob", HAPLOID, LIBRARY,  !ROOT, PAIR, 0, s+g+l, -1, 0.0f, "SpermBob"},
          };

        test(graph, expected);        
    }

    BOOST_TEST_CONTEXT("graph=clone_graph_autosomal") {
        libraries_t libs = {
            {"Dad", "Mom", "Bob", "Eve"},
            {"Dad", "Mom", "Bob", "Eve"}
        };

        Pedigree ped;
        ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"Dad"}});
        ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"Mom"}});
        ped.AddMember({"Bob",{"clone"},std::string{"Dad"},{},{},{},Sex::Male,{"Bob"}});
        ped.AddMember({"Eve",{"clone"},{},{},std::string{"Mom"},{},Sex::Female,{"Eve"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, false));

        const expected_graph_nodes_t expected = {
            {"GL/Dad",      DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/Mom",      DIPLOID, GERMLINE,  ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/Dad",      DIPLOID, LIBRARY,  !ROOT, PAIR, 0, s+l, -1, 0.0f, "Dad"},
            {"LB/Mom",      DIPLOID, LIBRARY,  !ROOT, PAIR, 1, s+l, -1, 0.0f, "Mom"},
            {"LB/Bob",      DIPLOID, LIBRARY,  !ROOT, PAIR, 0, 2*s+l, -1, 0.0f, "Bob"},
            {"LB/Eve",      DIPLOID, LIBRARY,  !ROOT, PAIR, 1, 2*s+l, -1, 0.0f, "Eve"},
          };

        test(graph, expected);        
    }
    BOOST_TEST_CONTEXT("graph=ploidy_graph") {
        libraries_t libs = {
            {"A", "B", "C", "D", "E", "F", "G"},
            {"A", "B", "C", "D", "E", "F", "G"},
        };

        Pedigree ped;
        ped.AddMember({"A",{"haploid"},  {},{},{},{},Sex::Autosomal,{"A"}});
        ped.AddMember({"B",{"gamete"},   {},{},{},{},Sex::Autosomal,{"B"}});
        ped.AddMember({"C",{"p=1"},      {},{},{},{},Sex::Autosomal,{"C"}});
        ped.AddMember({"D",{"ploidy=1"}, {},{},{},{},Sex::Autosomal,{"D"}});
        ped.AddMember({"E",{"diploid"},  {},{},{},{},Sex::Autosomal,{"E"}});
        ped.AddMember({"F",{"p=2"},      {},{},{},{},Sex::Autosomal,{"F"}});
        ped.AddMember({"G",{"ploidy=2"}, {},{},{},{},Sex::Autosomal,{"G"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_graph_nodes_t expected = {
            {"GL/A", HAPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/B", HAPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/C", HAPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/D", HAPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/E", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/F", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"GL/G", DIPLOID, GERMLINE, ROOT, FOUNDER, -1, 0.0f, -1, 0.0f, ""},
            {"LB/A", HAPLOID, LIBRARY, !ROOT, PAIR,     0, s+l,  -1, 0.0f, "A"},
            {"LB/B", HAPLOID, LIBRARY, !ROOT, PAIR,     1, s+l,  -1, 0.0f, "B"},
            {"LB/C", HAPLOID, LIBRARY, !ROOT, PAIR,     2, s+l,  -1, 0.0f, "C"},
            {"LB/D", HAPLOID, LIBRARY, !ROOT, PAIR,     3, s+l,  -1, 0.0f, "D"},
            {"LB/E", DIPLOID, LIBRARY, !ROOT, PAIR,     4, s+l,  -1, 0.0f, "E"},
            {"LB/F", DIPLOID, LIBRARY, !ROOT, PAIR,     5, s+l,  -1, 0.0f, "F"},
            {"LB/G", DIPLOID, LIBRARY, !ROOT, PAIR,     6, s+l,  -1, 0.0f, "G"},
          };

        test(graph, expected);        
    }
}

BOOST_AUTO_TEST_CASE(test_RelationshipGraph_Construct_exceptions) {
    using Sex = Pedigree::Sex;
    libraries_t libs = {
        {"A", "B", "C", "D", "E", "F", "G"},
        {"A", "B", "C", "D", "E", "F", "G"},
    };

    RelationshipGraph graph;
    std::string A = "A", B = "B", C = "C", D = "D", E="E";
    constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

    auto test = [&](Pedigree& pedigree, InheritanceModel model) {
        graph.Construct(pedigree, libs, model, g, s, l, true);
    };

    Pedigree clone_two_parents{{{A,{"clone"},B,{},C,{},Sex::Male,{A}}}};
    BOOST_CHECK_THROW(test(clone_two_parents,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree clone_bad_parent{{
        {A,{},{},{},{},{},Sex::Male,{A}},
        {B,{"clone"},C,{},{},{},Sex::Male,{B}}
    }};
    BOOST_CHECK_THROW(test(clone_bad_parent,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree gamete_two_parents{{{A,{"gamete"},B,{},C,{},Sex::Male,{A}}}};
    BOOST_CHECK_THROW(test(gamete_two_parents,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree gamete_bad_parent{{
        {A,{},{},{},{},{},Sex::Male,{A}},
        {B,{"gamete"},C,{},{},{},Sex::Male,{B}}
    }};
    BOOST_CHECK_THROW(test(gamete_bad_parent,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree gamete_bad_dad{{
        {A,{},{},{},{},{},Sex::Female,{A}},
        {B,{"gamete"},A,{},{},{},Sex::Male,{B}}
    }};
    BOOST_CHECK_THROW(test(gamete_bad_dad,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree gamete_bad_mom{{
        {A,{},{},{},{},{},Sex::Male,{A}},
        {B,{"gamete"},{},{},A,{},Sex::Male,{B}}
    }};
    BOOST_CHECK_THROW(test(gamete_bad_mom,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree no_dad{{
        {A,{},{},{},{},{},Sex::Male,  {A}},
        {B,{},{},{},{},{},Sex::Female,{B}},
        {C,{},{},{},B,{},Sex::Male,  {C}}
    }};
    BOOST_CHECK_THROW(test(no_dad,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree no_mom{{
        {A,{},{},{},{},{},Sex::Male,  {A}},
        {B,{},{},{},{},{},Sex::Female,{B}},
        {C,{},A,{},{},{},Sex::Male,  {C}}
    }};
    BOOST_CHECK_THROW(test(no_mom,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree bad_dad{{
        {A,{},{},{},{},{},Sex::Male,  {A}},
        {B,{},{},{},{},{},Sex::Female,{B}},
        {C,{},D,{},B,{},Sex::Male,  {C}}
    }};
    BOOST_CHECK_THROW(test(bad_dad,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree bad_mom{{
        {A,{},{},{},{},{},Sex::Male,  {A}},
        {B,{},{},{},{},{},Sex::Female,{B}},
        {C,{},A,{},D,{},Sex::Male,  {C}}
    }};
    BOOST_CHECK_THROW(test(bad_mom,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree female_dad{{
        {A,{},{},{},{},{},Sex::Female,  {A}},
        {B,{},{},{},{},{},Sex::Female,{B}},
        {C,{},A,{},B,{},Sex::Male,  {C}}
    }};
    BOOST_CHECK_THROW(test(female_dad,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree male_mom{{
        {A,{},{},{},{},{},Sex::Male,  {A}},
        {B,{},{},{},{},{},Sex::Male,{B}},
        {C,{},A,{},B,{},Sex::Male,  {C}}
    }};
    BOOST_CHECK_THROW(test(male_mom,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree selfing{{
        {A,{},{},{},{},{},Sex::Autosomal,  {A}},
        {B,{},{},{},{},{},Sex::Autosomal,{B}},
        {C,{},A,{},A,{},Sex::Autosomal,  {C}}
    }};
    BOOST_CHECK_THROW(test(selfing,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree empty_soma{{
        {A,{},{},{},{},{},Sex::Autosomal, {""}},
    }};
    BOOST_CHECK_THROW(test(empty_soma,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree bad_soma{{
        {A,{},{},{},{},{},Sex::Autosomal, {"(A;"}},
    }};
    BOOST_CHECK_THROW(test(bad_soma,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree loop{{
        {A,{},{},{},{},{},Sex::Autosomal,{A}},
        {B,{},{},{},{},{},Sex::Autosomal,{B}},
        {C,{},A,{},B,{},  Sex::Autosomal,{C}},
        {D,{},A,{},B,{},  Sex::Autosomal,{D}},
        {E,{},C,{},D,{},  Sex::Autosomal,{E}}
    }};
    BOOST_CHECK_THROW(test(loop,InheritanceModel::Autosomal), std::invalid_argument);

    Pedigree bad_sexlinkage{{
        {A,{},{},{},{},{},Sex::Autosomal,{A}},
        {B,{},{},{},{},{},Sex::Autosomal,{B}},
        {C,{},A,{},B,{},Sex::Male,{A}},
        {D,{},B,{},A,{},Sex::Female,{B}},
    }};
    BOOST_CHECK_THROW(test(bad_sexlinkage,InheritanceModel::XLinked), std::invalid_argument);
    BOOST_CHECK_THROW(test(bad_sexlinkage,InheritanceModel::YLinked), std::invalid_argument);
    BOOST_CHECK_THROW(test(bad_sexlinkage,InheritanceModel::ZLinked), std::invalid_argument);
    BOOST_CHECK_THROW(test(bad_sexlinkage,InheritanceModel::WLinked), std::invalid_argument);
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
    quad_ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"DadSm"}});
    quad_ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"MomSm"}});
    quad_ped.AddMember({"Eve",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Female,{"EveSm"}});
    quad_ped.AddMember({"Bob",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Male,{"BobSm"}});

    BOOST_TEST_CONTEXT("graph=quad_graph_autosomal") {
        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {1,3}},
            {Op::TOFATHERFAST, {0,1,4,5}},
            {Op::UP, {0,2}},
        };

        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::Autosomal, g, s, l, true));
        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_xlinked") {
        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {1,3}},
            {Op::UP, {1,5}},
            {Op::TOFATHERFAST, {0,1,4}},
            {Op::UP, {0,2}},
        };

        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::XLinked, g, s, l, true));
        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_ylinked") {
        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {0,1}},
            {Op::UP, {0,2}},
        };

        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::YLinked, g, s, l, true));
        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_zlinked") {
        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {1,3}},
            {Op::TOFATHERFAST, {0,1,5}},
            {Op::UP, {0,2}},
            {Op::UP, {0,4}},
        };

        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::ZLinked, g, s, l, true));
        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_wlinked") {
        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {0,1}},
            {Op::UP, {0,2}},
        };

        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::WLinked, g, s, l, true));
        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_maternal") {
        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {0,2}},
            {Op::UPFAST, {1,3}},
            {Op::UP, {1,5}},
            {Op::UP, {1,4}},
        };

        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::Maternal, g, s, l, true));
        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=quad_graph_paternal") {
        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {0,2}},
            {Op::UP, {0,5}},
            {Op::UP, {0,4}},
            {Op::UPFAST, {1,3}},
        };

        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(quad_ped, quad_libs,
            InheritanceModel::Paternal, g, s, l, true));
        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=somatic_graph_normalized") {
        libraries_t somatic_libs = {
            {"M1A", "M1B", "M2", "M3A", "M3B"},
            {"M1", "M1", "M2", "M3", "M3"}
        };

        Pedigree somatic_ped;
        somatic_ped.AddMember({"M",{},{},{},{},{},Sex::Female,{"((M1,M2)M12,M3);"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(somatic_ped, somatic_libs,
            InheritanceModel::Autosomal, g, s, l, true));

        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {2,8}},
            {Op::UP,     {2,9}},
            {Op::UPFAST, {1,2}},
            {Op::UPFAST, {4,5}},
            {Op::UP,     {4,6}},
            {Op::UPFAST, {3,4}},
            {Op::UP,     {3,7}},
            {Op::UP,     {1,3}},
            {Op::UPFAST, {0,1}},
        };

        test(graph, expected);
    }

   BOOST_TEST_CONTEXT("graph=twofam_graph") {
        libraries_t twofam_libs = {
            {"DadLb", "MomLb","EveLb", "BobLb", "SamLb"},
            {"DadSm", "MomSm", "EveSm", "BobSm", "SamSm"}
        };

        Pedigree twofam_ped;
        twofam_ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"DadSm"}});
        twofam_ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"MomSm"}});
        twofam_ped.AddMember({"Eve",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Female,{"EveSm"}});
        twofam_ped.AddMember({"Bob",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Male,{"BobSm"}});
        twofam_ped.AddMember({"Sam",{},{},{},{},{},Sex::Male,{"SamSm"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(twofam_ped, twofam_libs,
            InheritanceModel::Autosomal, g, s, l, true));

        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {1,4}},
            {Op::TOFATHERFAST, {0,1,5,6}},
            {Op::UP, {0,3}},
            {Op::UPFAST, {2,7}}
        };

        test(graph, expected);
    }

    BOOST_TEST_CONTEXT("graph=m12") {
/*
1-2    3-4
 |      |
 7------8    5-6
   |  |       |
   9  10-----11
          |
          12
*/
        libraries_t libs = {
            {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"},
            {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"}
        };

        Pedigree ped;
        ped.AddMember({"1",{},{},{},{},{},Sex::Male,{"1"}});
        ped.AddMember({"2",{},{},{},{},{},Sex::Female,{"2"}});
        ped.AddMember({"3",{},{},{},{},{},Sex::Male,{"3"}});
        ped.AddMember({"4",{},{},{},{},{},Sex::Female,{"4"}});
        ped.AddMember({"5",{},{},{},{},{},Sex::Male,{"5"}});
        ped.AddMember({"6",{},{},{},{},{},Sex::Female,{"6"}});
        ped.AddMember({"7",{},std::string{"1"},{},std::string{"2"},{},Sex::Male,{"7"}});
        ped.AddMember({"8",{},std::string{"3"},{},std::string{"4"},{},Sex::Female,{"8"}});
        ped.AddMember({"9",{},std::string{"7"},{},std::string{"8"},{},Sex::Female,{"9"}});
        ped.AddMember({"10",{},std::string{"7"},{},std::string{"8"},{},Sex::Female,{"10"}});
        ped.AddMember({"11",{},std::string{"5"},{},std::string{"6"},{},Sex::Male,{"11"}});
        ped.AddMember({"12",{},std::string{"11"},{},std::string{"10"},{},Sex::Male,{"12"}});

        constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

        //Construct Graph
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(ped, libs, InheritanceModel::Autosomal, g, s, l, true));

        const expected_peeling_nodes_t expected = {
            {Op::UPFAST, {1,11}},
            {Op::UPFAST, {2,12}},
            {Op::UPFAST, {3,13}},
            {Op::TOCHILDFAST, {2,3,7}},
            {Op::UPFAST, {4,14}},
            {Op::UPFAST, {5,15}},
            {Op::TOCHILDFAST, {4,5,8}},       
            {Op::UPFAST, {8,20}},   
            {Op::TOMOTHERFAST, {8,9,21}},
            {Op::UP, {9,19}}, 
            {Op::UPFAST, {7,17}},
            {Op::TOFATHERFAST, {6,7,9,18}},
            {Op::UP, {6,16}},         
            {Op::TOFATHERFAST, {0,1,6}},         
            {Op::UP, {0,10}},         
        };

        test(graph, expected);
    }

}

BOOST_AUTO_TEST_CASE(test_RelationshipGraph_Peel) {
    using boost::copy;
    using boost::generate;
    using boost::sort;

    const double prec = 2*DBL_EPSILON;

    xorshift64 xrand(++g_seed_counter);

    using Op = peel::Op;

    constexpr float g = 1e-8, s = 3e-8, l = 4e-8;

    libraries_t twofam_libs = {
        {"DadLb", "MomLb","EveLb", "BobLb", "SamLb"},
        {"DadSm", "MomSm", "EveSm", "BobSm", "SamSm"}
    };

    Pedigree twofam_ped;
    twofam_ped.AddMember({"Dad",{},{},{},{},{},Sex::Male,{"DadSm"}});
    twofam_ped.AddMember({"Mom",{},{},{},{},{},Sex::Female,{"MomSm"}});
    twofam_ped.AddMember({"Eve",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Female,{"EveSm"}});
    twofam_ped.AddMember({"Bob",{},std::string{"Dad"},{},std::string{"Mom"},{},Sex::Male,{"BobSm"}});
    twofam_ped.AddMember({"Sam",{},{},{},{},{},Sex::Male,{"SamSm"}});


    BOOST_TEST_CONTEXT("graph=twofam_graph_autosomal") {
        RelationshipGraph graph;
        BOOST_REQUIRE_NO_THROW(graph.Construct(twofam_ped, twofam_libs,
            InheritanceModel::Autosomal, g, s, l, true));

        auto dmod = mutation::Model{1e-6, 4};
        auto mmod = mutation::Model{0.7e-6, 4};
        auto da = dmod.TransitionMatrix(4);
        auto ma = mmod.TransitionMatrix(4);

        auto work = graph.CreateWorkspace();

        // Create random log likelihoods sorted in decreasing order
        std::vector<std::vector<double>> expected_lib_lower(5);
        for(auto &&a : expected_lib_lower) {
            a.resize(10);
            generate(a, [&](){ return xrand.get_double52(); });
            sort(a, std::greater<double>());
        }
        for(int i = work.library_nodes.first, u=0;i<work.library_nodes.second;++i) {
            work.lower[i].resize(10);
            copy(expected_lib_lower[u++], work.lower[i].data());
        }
        // Create random priors sorted in decreasing order
        std::vector<std::vector<double>> expected_founder_upper(3);
        for(auto &&a : expected_founder_upper) {
            a.resize(10);
            generate(a, [&](){ return xrand.get_double52(); });
            sort(a, std::greater<double>());
        }
        for(int i = work.founder_nodes.first,u=0;i<work.founder_nodes.second;++i) {
            work.upper[i].resize(10);
            copy(expected_founder_upper[u++], work.upper[i].data());
        }
        // Create transition matrices
        TransitionMatrixVector mats(8);
        mats[3] = mutation::mitosis_matrix(4, dmod, mutation::transition_t{}, 2);
        mats[4] = mutation::mitosis_matrix(4, mmod, mutation::transition_t{}, 2);
        mats[5] = mutation::meiosis_matrix(4, dmod, mmod, mutation::transition_t{}, 2, 2);
        mats[6] = mutation::meiosis_matrix(4, dmod, mmod, mutation::transition_t{}, 2, 2);
        mats[3] = mutation::mitosis_matrix(4, dmod, mutation::transition_t{}, 2);
        mats[7] = mutation::mitosis_matrix(4, dmod, mutation::transition_t{}, 2);

        double test_value = graph.PeelForwards(work, mats);
        // Family 1
        double d1 = 0.0;
        std::vector<double> lower0(10,0.0);
        std::vector<double> lower1(10,0.0);
        for(int i=0;i<10;++i) {
            for(int j=0;j<10;++j) {
                lower0[i] += expected_lib_lower[0][j]*mats[3](i,j);
                lower1[i] += expected_lib_lower[1][j]*mats[4](i,j);
            }
        }
        for(int i=0;i<10;++i) {
            for(int j=0;j<10;++j) {
                int ij = i*10+j;
                double a = 0.0;
                for(int k=0;k<10;++k) {
                    a += expected_lib_lower[2][k]*mats[5](ij,k);
                }
                double b = 0.0;
                for(int k=0;k<10;++k) {
                    b += expected_lib_lower[3][k]*mats[6](ij,k);
                }
                d1 += a*b*(lower1[j]*expected_founder_upper[1][j])
                         *(lower0[i]*expected_founder_upper[0][i]);
            }
        }        

        // Family 2
        double d2 = 0.0;
        for(int i=0;i<10;++i) {
            for(int j=0;j<10;++j) {
                d2 += expected_lib_lower[4][j]*mats[7](i,j)*expected_founder_upper[2][i];
            }
        }
        double expected_value = log(d1)+log(d2);

        BOOST_CHECK_CLOSE_FRACTION(test_value, expected_value, prec);

        // Check that backwards peeling produces proper marginals
        graph.PeelBackwards(work, mats);

        std::vector<double> test_marginal;
        for(int n=0;n<work.lower.size();++n) {
            double d = (work.lower[n]*work.upper[n]).sum();
            test_marginal.push_back(d);
        }
        std::vector<double> test_marginal_super;
        for(int n=3;n<work.lower.size();++n) {
            double d = (work.lower[n]*(mats[n].transpose()*work.super[n].matrix()).array()).sum();
            test_marginal_super.push_back(d);
        }
        std::vector<double> expected_marginal(work.lower.size(), 1.0);
        std::vector<double> expected_marginal_super(work.lower.size()-3, 1.0);

        CHECK_CLOSE_RANGES(test_marginal, expected_marginal, prec);
        CHECK_CLOSE_RANGES(test_marginal_super, expected_marginal_super, prec);

        // Testing Forward Peeling with Dirty Lower.
        double test_value_2 = graph.PeelForwards(work, mats);
        BOOST_CHECK_CLOSE_FRACTION(test_value_2, expected_value, prec);
    }
}
