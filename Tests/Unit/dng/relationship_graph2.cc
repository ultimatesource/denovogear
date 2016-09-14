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


#define BOOST_TEST_MODULE dng::relationship_graph2

#include <dng/relationship_graph.h>

#include <iostream>

#include "../boost_test_helper.h"
#include "fixture_read_trio_from_file.h"
#include "relationship_graph_helper.h"

struct FixturePedigreeMid {

    std::string fixture;

    dng::io::Pedigree io_pedigree;
    dng::ReadGroups rgs;

    typedef dng::task::Call task_type;
    struct arg_t : public task_type::argument_type {
        bool help;
        bool version;
        std::string arg_file;

        std::string run_name;
        std::string run_path;
    } arg;

    FixturePedigreeMid(std::string s = "FixturePedigreeMid") : fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        po::options_description ext_desc, int_desc;
        po::positional_options_description pos_desc;
        po::variables_map vm;

        int argc=4;
        char *argv[argc];
        argv[0] = (char*) "test";
        argv[1] = (char*) "-p";

        std::string ped_filename (TESTDATA_DIR);
        ped_filename.append("/relationship_graph/relationship_graph.ped");
        argv[2] = (char*) ped_filename.data();

        std::string vcf_filename = TESTDATA_DIR;
        vcf_filename.append("/relationship_graph/relationship_graph.vcf");
        argv[3] = (char*) vcf_filename.data();

        add_app_args(ext_desc,
                     static_cast<typename task_type::argument_type &>(arg));
        int_desc.add_options()("input",
                               po::value<std::vector<std::string> >(&arg.input),
                               "input files");
        int_desc.add(ext_desc);
        pos_desc.add("input", -1);
        po::store(
                po::command_line_parser(argc, argv).options(int_desc).positional(
                        pos_desc).run(), vm);
        po::notify(vm);

        // Parse pedigree from file
        std::ifstream ped_file(arg.ped);
        io_pedigree.Parse(istreambuf_range(ped_file));

        std::vector<hts::File> indata;
        std::vector<hts::bcf::File> bcfdata;
        for (auto &&str : arg.input) {
            indata.emplace_back(str.c_str(), "r");
            if (indata.back().is_open()) {
                continue;
            }
            throw std::runtime_error("unable to open input file '" + str + "'.");
        }
        bcfdata.emplace_back(std::move(indata[0]));
        rgs.ParseSamples(bcfdata[0]);
        arg.mu= 0.05;
        arg.mu_somatic = 0.07;
        arg.mu_library = 0.11;

    }

    ~FixturePedigreeMid() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }

};

/*
1-2    3-4
 |      |
 7------8   5-6
   |  |      |
   10 11-----9
          |
          12
*/

namespace dng {
BOOST_FIXTURE_TEST_CASE(test_constructor, FixturePedigreeMid ) {

    dng::RelationshipGraph relationship_graph;
    relationship_graph.Construct(io_pedigree, rgs, arg.mu, arg.mu_somatic, arg.mu_library);

    BOOST_CHECK_EQUAL(22, relationship_graph.num_nodes() );

    auto workspace = relationship_graph.CreateWorkspace();
    BOOST_CHECK_EQUAL(0, workspace.founder_nodes.first);
    BOOST_CHECK_EQUAL(6, workspace.founder_nodes.second);
    BOOST_CHECK_EQUAL(0, workspace.germline_nodes.first);
    BOOST_CHECK_EQUAL(10, workspace.germline_nodes.second);
    BOOST_CHECK_EQUAL(10, workspace.somatic_nodes.first);
    BOOST_CHECK_EQUAL(10, workspace.somatic_nodes.second);
    BOOST_CHECK_EQUAL(10, workspace.library_nodes.first);
    BOOST_CHECK_EQUAL(22, workspace.library_nodes.second);

    auto labels = relationship_graph.labels();

    const std::vector<std::string> expected_labels = {
            "GL-1", "GL-2",
            "GL-4", "GL-5",
            "GL-9", "GL-10",
            "GL-3", "GL-6",
            "GL-11", "GL-8",
            "LB-NA12001:Solexa-001", "LB-NA12002:Solexa-002",
            "LB-NA12003:Solexa-003", "LB-NA12004:Solexa-004",
            "LB-NA12005:Solexa-005", "LB-NA12006:Solexa-006",
            "LB-NA12007:Solexa-007", "LB-NA12008:Solexa-008",
            "LB-NA12009:Solexa-009", "LB-NA12010:Solexa-010",
            "LB-NA12011:Solexa-011", "LB-NA12012:Solexa-012"
    };
    boost_check_equal_vector(expected_labels, labels);


    auto transitions = relationship_graph.transitions();
    auto s_max = static_cast<size_t>(-1);
    double mu_somatic_library = arg.mu + arg.mu_somatic + arg.mu_library;
    double somatic_library = arg.mu_somatic + arg.mu_library;
    std::vector<RelationshipGraph::transition_t> expected_transitions = {
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},

            {RelationshipGraph::TransitionType::Germline, 0, 1, arg.mu, arg.mu},
            {RelationshipGraph::TransitionType::Germline, 2, 3, arg.mu, arg.mu},
            {RelationshipGraph::TransitionType::Germline, 5, 4, arg.mu, arg.mu},
            {RelationshipGraph::TransitionType::Germline, 6, 7, arg.mu, arg.mu},

            {RelationshipGraph::TransitionType::Somatic, 0, s_max,
                    somatic_library, 0},
			{RelationshipGraph::TransitionType::Somatic, 1, s_max,
			        somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 6, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 2, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 3, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 7, s_max,
                    somatic_library, 0},

            {RelationshipGraph::TransitionType::Germline, 6, 7,
                    mu_somatic_library, mu_somatic_library},
            {RelationshipGraph::TransitionType::Somatic, 9, s_max,
                    somatic_library, 0},

            {RelationshipGraph::TransitionType::Somatic, 4, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 5, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 8, s_max,
                    somatic_library, 0},

            {RelationshipGraph::TransitionType::Germline, 8, 9,
                    mu_somatic_library, mu_somatic_library}
    };

    BOOST_CHECK_EQUAL(expected_transitions.size(), transitions.size());
    for (int k = 0; k < expected_transitions.size(); ++k) {
        auto expected = expected_transitions[k];
        auto actual = transitions[k];
        BOOST_CHECK(expected.type == actual.type);
        BOOST_CHECK_EQUAL(expected.parent1, actual.parent1);
        BOOST_CHECK_EQUAL(expected.parent2, actual.parent2);
        BOOST_CHECK_CLOSE(expected.length1, actual.length1, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
        BOOST_CHECK_CLOSE(expected.length2, actual.length2, BOOST_CLOSE_PERCENTAGE_THRESHOLD);

    }

}


BOOST_FIXTURE_TEST_CASE(test_pedigree_inspect, FixturePedigreeMid) {

    dng::RelationshipGraph relationship_graph;
    relationship_graph.Construct(io_pedigree, rgs, arg.mu, arg.mu_somatic, arg.mu_library);
    BOOST_CHECK_EQUAL(22, relationship_graph.num_nodes());

    std::vector<peel::family_members_t> family = relationship_graph.family_members_;
    std::vector<peel::family_members_t> expected_family = {
            {2, 13},
            {3, 14},
            {2, 3, 7},
            {5, 19},
            {4, 18},
            {5, 4, 8},
            {8, 20},
            {8, 9, 21},
            {9, 17},
            {7, 15},
            {6, 7, 9, 16},
            {6, 12},
            {1, 11},
            {0, 1, 6},
            {0, 10}
    };

    for (int f = 0; f < expected_family.size(); ++f) {
        boost_check_equal_vector(expected_family[f], family[f]);
    }

	std::vector<decltype(peel::op::NUM)> ops = relationship_graph.peeling_ops_;
    std::vector<decltype(peel::op::NUM)> expected_ops = {
            peel::op::UP,
            peel::op::UP,
			peel::op::TOCHILD,
			peel::op::UP,
			peel::op::UP,
			peel::op::TOCHILD,
			peel::op::UP,
            peel::op::TOMOTHER,
            peel::op::UP,
            peel::op::UP,
            peel::op::TOFATHER,
            peel::op::UP,
            peel::op::UP,
            peel::op::TOFATHER,
            peel::op::UP
    };
	boost_check_equal_vector(expected_ops, ops);

	std::vector<decltype(peel::op::NUM)> functions_ops =
			relationship_graph.peeling_functions_ops_;
    std::vector<decltype(peel::op::NUM)> expected_functions_ops = {
            peel::op::UPFAST,
            peel::op::UPFAST,
            peel::op::TOCHILDFAST,
            peel::op::UPFAST,
            peel::op::UPFAST,
            peel::op::TOCHILDFAST,
            peel::op::UPFAST,
            peel::op::TOMOTHERFAST,
            peel::op::UP,
            peel::op::UPFAST,
            peel::op::TOFATHERFAST,
            peel::op::UP,
            peel::op::UPFAST,
            peel::op::TOFATHERFAST,
            peel::op::UP
    };
    boost_check_equal_vector(expected_functions_ops, functions_ops);

}


} // namespace dng





