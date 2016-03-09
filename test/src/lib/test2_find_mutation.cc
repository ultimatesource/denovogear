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


#define BOOST_TEST_MODULE dng::lib::find_mutation

#include <ctime>

#include <iostream>
#include <fstream>

#include <boost/test/unit_test.hpp>

#include <dng/task/call.h>
#include <dng/find_mutation.h>
#include <boost_test_helper.h>
#include <fixture/fixture_trio_workspace.h>



namespace utf = boost::unit_test;

const int NUM_TEST = 100;


// TODO: Example of BOOST_DATA_TEST_CASE and BOOST_PARAM_TEST_CASE.
// TODO: Should be able to replace the for loop with these.
// TODO: Might not be able to use fixture.
//
//
//std::vector<int> test_types;//
//
//BOOST_AUTO_TEST_SUITE(suite1,
//  * utf::fixture<Fx>(std::string("FX"))
//  * utf::fixture<Fx>(std::string("FX2")))
//
//  BOOST_AUTO_TEST_CASE(test1, * utf::fixture(&setup, &teardown){
//    BOOST_TEST(true);
//  }

//
//BOOST_DATA_TEST_CASE( test_case_arity1, data::xrange(5), my_var )
//{
//    BOOST_TEST_MESSAGE("running data: ");
//    BOOST_TEST((my_var <= 4 && my_var >= 0));
//}
//BOOST_PARAM_TEST_CASE(test_function, params_begin, params_end);
//BOOST_AUTO_TEST_SUITE_END()

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }


//BOOST_AUTO_TEST_SUITE(test_peeling_suite,  * utf::fixture<Fx>(std::string("FX")) )
//BOOST_FIXTURE_TEST_SUITE(test_find_mutation_suite, TrioWorkspace )


//BOOST_AUTO_TEST_CASE(test_constructor, *utf::fixture(&setup, &teardown)) {
BOOST_FIXTURE_TEST_CASE(test_constructor, TrioWorkspace,
                        *utf::fixture(&setup, &teardown)) {

    FindMutations find_mutation{min_prob, pedigree, test_param_1};

    BOOST_CHECK_EQUAL(find_mutation.min_prob_, 0.01);
    BOOST_CHECK_EQUAL(find_mutation.params_.theta, 0.001);
    BOOST_CHECK_EQUAL(find_mutation.params_.ref_weight, 1);

    std::array<double, 4> expect_freqs{0.3, 0.2, 0.2, 0.3};
    auto freqs1 = find_mutation.params_.nuc_freq;
    BoostCheckCloseVector(expect_freqs, freqs1, 4);

    std::vector<std::array<double, 4>> expect_gamma{{0.98, 0.0005, 0.0005, 1.04},
                                                    {0.02, 0.075, 0.005, 1.18}};
    auto gamma_a = find_mutation.params_.params_a;
    BOOST_CHECK_EQUAL(gamma_a.pi, expect_gamma[0][0]);
    BOOST_CHECK_EQUAL(gamma_a.phi, expect_gamma[0][1]);
    BOOST_CHECK_EQUAL(gamma_a.epsilon, expect_gamma[0][2]);
    BOOST_CHECK_EQUAL(gamma_a.omega, expect_gamma[0][3]);

    auto gamma_b = find_mutation.params_.params_b;
    BOOST_CHECK_EQUAL(gamma_b.pi, expect_gamma[1][0]);
    BOOST_CHECK_EQUAL(gamma_b.phi, expect_gamma[1][1]);
    BOOST_CHECK_EQUAL(gamma_b.epsilon, expect_gamma[1][2]);
    BOOST_CHECK_EQUAL(gamma_b.omega, expect_gamma[1][3]);

#if CALCULATE_ENTROPY == 1
    auto event = find_mutation.event_;
    BOOST_CHECK_EQUAL(event.size(), 5);
    for (int k = 0; k < 5; ++k) {
        BOOST_CHECK_EQUAL(event[k], 0);
    }
#endif
}



BOOST_FIXTURE_TEST_CASE(test_prior, TrioWorkspace, *utf::fixture(&setup, &teardown)) {

    std::vector<std::array<double, 10>> expected_prior {
        {0.998951118846171, 0.000199760259730275, 0.000199760259730275, 0.000299640389595412, 9.98701448476561e-05, 3.99400699250774e-08, 5.99101048876161e-08, 9.98701448476561e-05, 5.99101048876161e-08, 0.000149820194797706},
        {0.000149820194797706, 0.000299610434542968, 5.99101048876161e-08, 8.98651573314242e-08, 0.998801318621409, 0.000199740289695312, 0.000299610434542968, 9.98701448476561e-05, 5.99101048876161e-08, 0.000149820194797706},
        {0.000149820194797706, 5.99101048876161e-08, 0.000299610434542968, 8.98651573314242e-08, 9.98701448476561e-05, 0.000199740289695312, 5.99101048876161e-08, 0.998801318621409, 0.000299610434542968, 0.000149820194797706},
        {0.000149820194797706, 5.99101048876161e-08, 5.99101048876161e-08, 0.000299640389595412, 9.98701448476561e-05, 3.99400699250774e-08, 0.000199760259730275, 9.98701448476561e-05, 0.000199760259730275, 0.998951118846171},
        {0.29979020979021, 0.00011988011988012, 0.00011988011988012, 0.00017982017982018, 0.19984015984016, 7.99200799200799e-05, 0.00011988011988012, 0.19984015984016, 0.00011988011988012, 0.29979020979021 }
    };

    FindMutations find_mutation {min_prob, pedigree, test_param_1};
    auto *pArray = find_mutation.genotype_prior_;

    for (int i = 0; i < 5; ++i) {
        auto prior_array = pArray[i];
        BoostCheckCloseVector(expected_prior[i], prior_array);
    }


}

BOOST_FIXTURE_TEST_CASE(test_full_transition, TrioWorkspace, *utf::fixture(&setup, &teardown)) {

    std::array<double, 4> freqs{0.3, 0.2, 0.2, 0.3};
    double mu = 1e-8;
    auto dad = f81::matrix(3*mu, freqs);
    auto mom = f81::matrix(3*mu, freqs);

    auto exp_germline_full = meiosis_diploid_matrix(dad, mom);
    auto exp_germline_nomut = meiosis_diploid_matrix(dad, mom, 0);
    auto exp_germline_posmut = exp_germline_full - exp_germline_nomut;
    auto exp_germline_onemut = meiosis_diploid_matrix(dad, mom, 1);
    auto exp_germline_mean = meiosis_diploid_mean_matrix(dad, mom);


    auto orig = f81::matrix(2*mu, freqs);
    auto exp_somatic_full = mitosis_diploid_matrix(orig);
    auto exp_somatic_nomut = mitosis_diploid_matrix(orig, 0);
    auto exp_somatic_posmut = exp_somatic_full - exp_somatic_nomut;
    auto exp_somatic_onemut = mitosis_diploid_matrix(orig, 1);
    auto exp_somatic_mean = mitosis_diploid_mean_matrix(orig);

    std::vector<TransitionMatrix> exp_full {{},{},exp_germline_full, exp_somatic_full, exp_somatic_full};
    std::vector<TransitionMatrix> exp_nomut {{},{},exp_germline_nomut, exp_somatic_nomut, exp_somatic_nomut};
    std::vector<TransitionMatrix> exp_posmut {{},{},exp_germline_posmut, exp_somatic_posmut, exp_somatic_posmut};
    std::vector<TransitionMatrix> exp_onemut {{},{},exp_germline_onemut, exp_somatic_onemut, exp_somatic_onemut};
    std::vector<TransitionMatrix> exp_mean {{},{},exp_germline_mean, exp_somatic_mean, exp_somatic_mean};

    FindMutations find_mutation {min_prob, pedigree, test_param_1};
    auto full_matrices = find_mutation.full_transition_matrices_;
    auto nomut_matrices = find_mutation.nomut_transition_matrices_;
    auto posmut_matrices = find_mutation.posmut_transition_matrices_;
    auto onemut_matrices = find_mutation.onemut_transition_matrices_;
    auto mean_matrices = find_mutation.mean_matrices_;
    for (int i = 0; i < 5; ++i) {
        BoostCheckMatrix(exp_full[i], full_matrices[i]);
        BoostCheckMatrix(exp_nomut[i], nomut_matrices[i]);
        BoostCheckMatrix(exp_posmut[i], posmut_matrices[i]);
        BoostCheckMatrix(exp_onemut[i], onemut_matrices[i]);
        BoostCheckMatrix(exp_mean[i], mean_matrices[i]);
    }
}


BOOST_FIXTURE_TEST_CASE(test_operator, TrioWorkspace, *utf::fixture(&setup, &teardown)) {


    std::vector<peel::family_members_t> family {
            {1,4},
            {0,1,2},
            {0,3}
    };
    MutationStats stats(min_prob);
    FindMutations find_mutation{min_prob, pedigree, test_param_1};
    find_mutation.CalculateMutation(read_depths, ref_index, stats);


    double scale = setup_workspace(ref_index, read_depths);

    //Test basic stats
    auto nomut_matrices = find_mutation.nomut_transition_matrices_;
    workspace.CleanupFast();
    dng::peel::up(workspace, family[0], nomut_matrices);
    dng::peel::to_father_fast(workspace, family[1], nomut_matrices);
    dng::peel::up(workspace, family[2], nomut_matrices);
    double result_nomut = log((workspace.lower[0] * workspace.upper[0]).sum());


    auto full_matrices = find_mutation.full_transition_matrices_;
    workspace.CleanupFast();
    dng::peel::up(workspace, family[0], full_matrices );
    dng::peel::to_father(workspace, family[1], full_matrices );
    dng::peel::up(workspace, family[2], full_matrices );
    double result_full = log((workspace.lower[0] * workspace.upper[0]).sum());


    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    const double expected_mup = -std::expm1(result_nomut - result_full);
    const double expected_lld = (result_full +scale) / M_LN10;
    const double expected_llh = result_full / M_LN10;

    BOOST_CHECK_CLOSE(expected_mup, stats.mup_, BOOST_CLOSE_THRESHOLD);
    BOOST_CHECK_CLOSE(expected_lld, stats.lld_, BOOST_CLOSE_THRESHOLD);
    BOOST_CHECK_CLOSE(expected_llh, stats.llh_, BOOST_CLOSE_THRESHOLD);

    //Test posterior
    full_matrices = find_mutation.full_transition_matrices_;
    workspace.CleanupFast();
    dng::peel::up(workspace, family[0], full_matrices );
    dng::peel::to_father(workspace, family[1], full_matrices );
    dng::peel::up(workspace, family[2], full_matrices );

    dng::peel::up_reverse(workspace, family[2], full_matrices );
    dng::peel::to_father_reverse(workspace, family[1], full_matrices );
    dng::peel::up_reverse(workspace, family[0], full_matrices );

    dng::GenotypeArray expected_posterior;
    for (std::size_t i = 0; i < workspace.num_nodes; ++i) {
        expected_posterior = WORKSPACE_T_MULTIPLE_UPPER_LOWER(workspace, i);
        expected_posterior /= expected_posterior.sum();

        BoostCheckCloseVector(expected_posterior, stats.posterior_probabilities_[i]);
    }



}


BOOST_FIXTURE_TEST_CASE(test_calculate_mutation_expected, TrioWorkspace, *utf::fixture(&setup, &teardown)) {

    FindMutations find_mutation{min_prob, pedigree, test_param_1};

    MutationStats mutation_stats(min_prob);
    find_mutation.CalculateMutation(read_depths, ref_index, mutation_stats);

    BOOST_CHECK_CLOSE(0.116189, mutation_stats.mup_, BOOST_CLOSE_THRESHOLD);
    BOOST_CHECK_CLOSE(-30.5967, mutation_stats.lld_, BOOST_CLOSE_THRESHOLD);
    BOOST_CHECK_CLOSE(-6.68014, mutation_stats.llh_, BOOST_CLOSE_THRESHOLD);
    BOOST_CHECK_CLOSE(0.116189, mutation_stats.mux_, BOOST_CLOSE_THRESHOLD);
    BOOST_CHECK_CLOSE(0.116189, mutation_stats.mu1p_, BOOST_CLOSE_THRESHOLD);

    BOOST_CHECK_EQUAL("GGxGG>GT", mutation_stats.dnt_);
    BOOST_CHECK_EQUAL("LB-NA12878:Solexa-135852", mutation_stats.dnl_);
    BOOST_CHECK_EQUAL(42, mutation_stats.dnq_);


}


BOOST_FIXTURE_TEST_CASE(test_calculate_mutation, TrioWorkspace, *utf::fixture(&setup, &teardown)) {
    //TODO: Compare operator() with CalculateMutation ONLY, NOT a real test


    min_prob = 0; //Pass everything
    FindMutations::stats_t stats;
    FindMutations find_mutation{min_prob, pedigree, test_param_1};
    find_mutation.old_operator(read_depths, ref_index, &stats);

    MutationStats mutation_stats(min_prob);
    find_mutation.CalculateMutation(read_depths, ref_index, mutation_stats);

    BOOST_CHECK_EQUAL(stats.mup, mutation_stats.mup_);
    BOOST_CHECK_EQUAL(stats.llh, mutation_stats.llh_);
    BOOST_CHECK_EQUAL(stats.lld, mutation_stats.lld_);
    BOOST_CHECK_EQUAL(stats.mux, mutation_stats.mux_);
    BOOST_CHECK_EQUAL(stats.has_single_mut, mutation_stats.has_single_mut_);
    BOOST_CHECK_EQUAL(stats.mu1p, mutation_stats.mu1p_);

    BOOST_CHECK_EQUAL(stats.dnt, mutation_stats.dnt_);
    BOOST_CHECK_EQUAL(stats.dnl, mutation_stats.dnl_);
    BOOST_CHECK_EQUAL(stats.dnq, mutation_stats.dnq_);

    for (int i = 0; i < stats.posterior_probabilities.size(); ++i) {
        BoostCheckCloseVector(stats.posterior_probabilities[i],
                              mutation_stats.posterior_probabilities_[i]);
    }
    for (int i = 2; i < stats.genotype_likelihoods.size(); ++i) {
        BoostCheckCloseVector(stats.genotype_likelihoods[i],
                              mutation_stats.genotype_likelihoods_[i]);
    }
    for (int i = 2; i < stats.node_mup.size(); ++i) {
        BOOST_CHECK_CLOSE(stats.node_mup[i], mutation_stats.node_mup_[i],
                          BOOST_CLOSE_THRESHOLD);
    }

    for (int i = 2; i < stats.node_mu1p.size(); ++i) {
        BOOST_CHECK_CLOSE(stats.node_mu1p[i], mutation_stats.node_mu1p_[i],
                          BOOST_CLOSE_THRESHOLD);
    }

#if CALCULATE_ENTROPY == 1
    BOOST_CHECK_EQUAL(stats.dnc, mutation_stats.dnc_);
#endif

}



//BOOST_FIXTURE_TEST_SUITE(test_find_mutation_suite, TrioWorkspace )
//BOOST_AUTO_TEST_CASE(test_find_mutation,  *utf::fixture(&setup, &teardown)) {}
//BOOST_AUTO_TEST_SUITE_END()
//With suite
//Test module "dng::lib::find_mutation" has passed with:
//6 test cases out of 6 passed
//6298 assertions out of 6298 passed
//
//        Test suite "test_find_mutation_suite" has passed with:
//6 test cases out of 6 passed
//6298 assertions out of 6298 passed
//
//        Test case "test_find_mutation_suite/test_constructor" has passed with:
//23 assertions out of 23 passed


//Without suite
//    Test module "dng::lib::find_mutation" has passed with:
//  6 test cases out of 6 passed
//  6298 assertions out of 6298 passed
//
//  Test case "test_constructor" has passed with:
//    23 assertions out of 23 passed

