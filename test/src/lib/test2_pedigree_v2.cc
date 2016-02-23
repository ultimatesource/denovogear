//
// Created by steven on 2/9/16.
//

#define BOOST_TEST_MODULE dng::lib::pedigree

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <dng/pedigree_v2.h>
#include <dng/pedigree.h>

#include <fixture/fixture_read_trio_from_file.h>
#include <boost_test_helper.h>


namespace utf = boost::unit_test;


struct FixturePedigree : public ReadTrioFromFile{


    dng::Pedigree pedigree;
    dng::PedigreeV2 pedigree_v2;

    FixturePedigree(std::string s = "FixturePedigree") : ReadTrioFromFile(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        pedigree.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library);
        pedigree_v2.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library);

    }

    ~FixturePedigree() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }

};

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }



BOOST_FIXTURE_TEST_SUITE(test_pedigree_suite, FixturePedigree )

/*

0     1
|     |
---|---
|  |  |
3  2  4

*/
BOOST_AUTO_TEST_CASE(test_constructor, *utf::fixture(&setup, &teardown)) {

    BOOST_CHECK_EQUAL(5, pedigree.num_nodes() );

    auto workspace = pedigree.CreateWorkspace();
    BOOST_CHECK_EQUAL(0, workspace.founder_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.founder_nodes.second);
    BOOST_CHECK_EQUAL(0, workspace.germline_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.germline_nodes.second);
    BOOST_CHECK_EQUAL(2, workspace.somatic_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.somatic_nodes.second);
    BOOST_CHECK_EQUAL(2, workspace.library_nodes.first);
    BOOST_CHECK_EQUAL(5, workspace.library_nodes.second);


    auto labels = pedigree.labels();

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

    auto transitions = pedigree.transitions();
    auto size_t_negative_one = static_cast<size_t>(-1);
    std::vector<Pedigree::transition_t> expected_transitions = {
            {Pedigree::TransitionType::Founder, size_t_negative_one, size_t_negative_one , 0, 0},
            {Pedigree::TransitionType::Founder, size_t_negative_one, size_t_negative_one , 0, 0},
            {Pedigree::TransitionType::Germline, 0, 1, arg.mu+arg.mu_somatic+arg.mu_library, arg.mu+arg.mu_somatic+arg.mu_library},
            {Pedigree::TransitionType::Somatic, 0, size_t_negative_one, arg.mu_somatic + arg.mu_library, 0},
            {Pedigree::TransitionType::Somatic, 1, size_t_negative_one, arg.mu_somatic + arg.mu_library,0}
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
        BOOST_CHECK_CLOSE(expected.length1, actual.length1, BOOST_CLOSE_THRESHOLD);
        BOOST_CHECK_CLOSE(expected.length2, actual.length2, BOOST_CLOSE_THRESHOLD);

    }




}



BOOST_AUTO_TEST_CASE(test_pedigree_equal, *utf::fixture(&setup, &teardown)) {

    bool is_equal = pedigree.Equal(pedigree_v2);
    BOOST_CHECK(is_equal);
}

BOOST_AUTO_TEST_CASE(test_constructor_2, *utf::fixture(&setup, &teardown)) {

    BOOST_CHECK_EQUAL(5, pedigree_v2.num_nodes() );

    auto workspace = pedigree_v2.CreateWorkspace();
    BOOST_CHECK_EQUAL(0, workspace.founder_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.founder_nodes.second);
    BOOST_CHECK_EQUAL(0, workspace.germline_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.germline_nodes.second);
    BOOST_CHECK_EQUAL(2, workspace.somatic_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.somatic_nodes.second);
    BOOST_CHECK_EQUAL(2, workspace.library_nodes.first);
    BOOST_CHECK_EQUAL(5, workspace.library_nodes.second);


    auto labels = pedigree_v2.labels();

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

    
    auto size_t_negative_one = static_cast<size_t>(-1);
    std::vector<Pedigree::transition_t> expected_transitions = {
            {Pedigree::TransitionType::Founder, size_t_negative_one, size_t_negative_one , 0, 0},
            {Pedigree::TransitionType::Founder, size_t_negative_one, size_t_negative_one , 0, 0},
            {Pedigree::TransitionType::Germline, 0, 1, arg.mu+arg.mu_somatic+arg.mu_library, arg.mu+arg.mu_somatic+arg.mu_library},
            {Pedigree::TransitionType::Somatic, 0, size_t_negative_one, arg.mu_somatic + arg.mu_library, 0},
            {Pedigree::TransitionType::Somatic, 1, size_t_negative_one, arg.mu_somatic + arg.mu_library,0}
    };

    auto transitions = pedigree_v2.transitions();
    for (int k = 0; k < 5; ++k) {
        auto expected = expected_transitions[k];
        auto actual = transitions[k];
        BOOST_CHECK(expected.type == actual.type);
        BOOST_CHECK_EQUAL(expected.parent1, actual.parent1);
        BOOST_CHECK_EQUAL(expected.parent2, actual.parent2);
        BOOST_CHECK_CLOSE(expected.length1, actual.length1, BOOST_CLOSE_THRESHOLD);
        BOOST_CHECK_CLOSE(expected.length2, actual.length2, BOOST_CLOSE_THRESHOLD);

    }

}




BOOST_AUTO_TEST_CASE(test_pedigree_v2, *utf::fixture(&setup, &teardown)) {

    
    BOOST_CHECK_EQUAL(5, pedigree_v2.num_nodes());


    std::vector<peel::family_members_t> family = pedigree_v2.inspect_family_members();
    std::vector<peel::family_members_t> expected_family = {
            {1, 4},
            {0, 1, 2},
            {0, 3}
    };

    for (int f = 0; f < expected_family.size(); ++f) {
        boost_check_equal_vector(expected_family, family);
    }

    std::vector<decltype(peel::op::NUM)> ops = pedigree_v2.inspect_peeling_ops();
    std::vector<decltype(peel::op::NUM)> expected_ops = {peel::op::UP, peel::op::TOFATHER, peel::op::UP};
    boost_check_equal_vector(expected_ops, ops);



    std::vector<decltype(peel::op::NUM)> functions_ops = pedigree_v2.inspect_peeling_functions_ops();
    std::vector<decltype(peel::op::NUM)> expected_functions_ops = {peel::op::UPFAST, peel::op::TOFATHERFAST,
                                                                  peel::op::UP};
    boost_check_equal_vector(expected_functions_ops, functions_ops);


}

BOOST_AUTO_TEST_SUITE_END()

