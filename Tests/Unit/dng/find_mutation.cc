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


#define BOOST_TEST_MODULE dng::find_mutation

#include <ctime>

#include <iostream>
#include <fstream>

#include <dng/task/call.h>
#include <dng/find_mutations.h>

#include "../boost_test_helper.h"
#include "fixture_trio_workspace.h"

namespace dng {

BOOST_FIXTURE_TEST_CASE(test_constructor, TrioWorkspace) {

    FindMutations find_mutation{min_prob, pedigree, test_param_1};

    BOOST_CHECK_EQUAL(find_mutation.min_prob_, 0.01);
    BOOST_CHECK_EQUAL(find_mutation.params_.theta, 0.001);
    BOOST_CHECK_EQUAL(find_mutation.params_.ref_weight, 1);

    std::array<double, 4> expect_freqs{0.3, 0.2, 0.2, 0.3};
    auto freqs1 = find_mutation.params_.nuc_freq;
    boost_check_close_vector(expect_freqs, freqs1);

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



BOOST_FIXTURE_TEST_CASE(test_prior, TrioWorkspace) {

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
        boost_check_close_vector(expected_prior[i], prior_array);
    }


}

BOOST_FIXTURE_TEST_CASE(test_full_transition, TrioWorkspace) {

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
        boost_check_matrix(exp_full[i], full_matrices[i]);
        boost_check_matrix(exp_nomut[i], nomut_matrices[i]);
        boost_check_matrix(exp_posmut[i], posmut_matrices[i]);
        boost_check_matrix(exp_onemut[i], onemut_matrices[i]);
        boost_check_matrix(exp_mean[i], mean_matrices[i]);
    }
}


BOOST_FIXTURE_TEST_CASE(test_operator, TrioWorkspace) {


    std::vector<peel::family_members_t> family {
            {1,4},
            {0,1,2},
            {0,3}
    };
    FindMutations::stats_t mutation_stats;
//    MutationStats stats(min_prob);
    FindMutations find_mutation{min_prob, pedigree, test_param_1};
    find_mutation(read_depths, ref_index, &mutation_stats);


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

    BOOST_CHECK_CLOSE(expected_mup, mutation_stats.mup, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
    BOOST_CHECK_CLOSE(expected_lld, mutation_stats.lld, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
//    BOOST_CHECK_CLOSE(expected_llh, mutation_stats.llh, BOOST_CLOSE_PERCENTAGE_THRESHOLD);

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
        expected_posterior = workspace.lower[i] * workspace.upper[i];
        expected_posterior /= expected_posterior.sum();

        boost_check_close_vector(expected_posterior,
        		mutation_stats.posterior_probabilities[i]);
    }



}


BOOST_FIXTURE_TEST_CASE(test_calculate_mutation_expected, TrioWorkspace) {

    FindMutations find_mutation{min_prob, pedigree, test_param_1};

    FindMutations::stats_t mutation_stats;
//    MutationStats mutation_stats(min_prob);
    find_mutation(read_depths, ref_index, &mutation_stats);

    BOOST_CHECK_SMALL(0.116189 - mutation_stats.mup, 1e-6);
    BOOST_CHECK_SMALL(-30.5967 - mutation_stats.lld, 1e-4);
//    BOOST_CHECK_SMALL(-6.68014 - mutation_stats.llh, 1e-5);
    BOOST_CHECK_SMALL(0.116189 - mutation_stats.mux, 1e-6);
    BOOST_CHECK_SMALL(0.116189 - mutation_stats.mu1p, 1e-6);

    BOOST_CHECK_EQUAL("GGxGG>GT", mutation_stats.dnt);
    BOOST_CHECK_EQUAL("LB-NA12878:Solexa-135852", mutation_stats.dnl);
    BOOST_CHECK_EQUAL(42, mutation_stats.dnq);

#if CALCULATE_ENTROPY == 1
    BOOST_CHECK_EQUAL(100, mutation_stats.dnc);
#endif

}

} // namespace dng
