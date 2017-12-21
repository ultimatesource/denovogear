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


#define BOOST_TEST_MODULE dng::log_probability

#include <dng/probability.h>

#include <iostream>
#include <fstream>
#include <numeric>

#include <dng/task/call.h>

#include "../testing.h"

namespace dng {
    struct unittest_dng_log_probability {
        GETTER2(LogProbability, params, theta)
        GETTER2(LogProbability, params, asc_bias_hom)
        GETTER2(LogProbability, params, asc_bias_het)
        GETTER2(LogProbability, params, asc_bias_hap)
        GETTER2(LogProbability, params, over_dispersion_hom)
        GETTER2(LogProbability, params, over_dispersion_het)
        GETTER2(LogProbability, params, ref_bias)
        GETTER2(LogProbability, params, error_rate)
        GETTER2(LogProbability, params, error_entropy)
        GETTER2(LogProbability, params, mutation_entropy)
    };
}

using u = dng::unittest_dng_log_probability;

using namespace dng;
using namespace dng::detail;
using Sex = dng::Pedigree::Sex;

// Use a lambda function to construct a global relationship graph
RelationshipGraph rel_graph = []() -> RelationshipGraph {
    libraries_t libs = {
        {"Mom", "Dad", "Eve"},
        {"Mom", "Dad", "Eve"}
    };
    Pedigree ped;
    ped.AddMember("Dad","0","0",Sex::Male,"");
    ped.AddMember("Mom","0","0",Sex::Female,"");
    ped.AddMember("Eve","Dad","Mom",Sex::Female,"");

    RelationshipGraph g;
    g.Construct(ped, libs, 1e-8, 1e-8, 1e-8,true);
    return g;
}();

BOOST_AUTO_TEST_CASE(test_constructor_1) {
    LogProbability::params_t params;
    params.theta = 0.001;
    params.asc_bias_hom = 0.0;
    params.asc_bias_het = 0.0;
    params.asc_bias_hap = 0.0;
    params.over_dispersion_hom = 1e-4;
    params.over_dispersion_het = 1e-4;
    params.ref_bias = 1;
    params.error_rate = 1e-4;
    params.error_entropy = log(4);
    params.mutation_entropy = log(4);

    LogProbability log_probability{rel_graph, params};

    BOOST_CHECK_EQUAL(u::theta(log_probability), params.theta);
    BOOST_CHECK_EQUAL(u::asc_bias_hom(log_probability), params.asc_bias_hom);
    BOOST_CHECK_EQUAL(u::asc_bias_het(log_probability), params.asc_bias_het);
    BOOST_CHECK_EQUAL(u::asc_bias_hap(log_probability), params.asc_bias_hap);
    
    BOOST_CHECK_EQUAL(u::over_dispersion_hom(log_probability), params.over_dispersion_hom);
    BOOST_CHECK_EQUAL(u::over_dispersion_het(log_probability), params.over_dispersion_het);
    BOOST_CHECK_EQUAL(u::ref_bias(log_probability), params.ref_bias);
    BOOST_CHECK_EQUAL(u::error_rate(log_probability), params.error_rate);
    BOOST_CHECK_EQUAL(u::error_entropy(log_probability), params.error_entropy);
    BOOST_CHECK_EQUAL(u::mutation_entropy(log_probability), params.mutation_entropy);
}

// BOOST_FIXTURE_TEST_CASE(test_full_transition, TrioWorkspace) {

//     // std::array<double, 4> freqs{0.3, 0.2, 0.2, 0.3};
//     // double mu = 1e-8;
//     // auto dad = f81::matrix(3*mu, freqs);
//     // auto mom = f81::matrix(3*mu, freqs);

//     // auto exp_germline_full = meiosis_diploid_matrix(dad, mom);
//     // auto exp_germline_nomut = meiosis_diploid_matrix(dad, mom, 0);
//     // auto exp_germline_posmut = exp_germline_full - exp_germline_nomut;
//     // auto exp_germline_onemut = meiosis_diploid_matrix(dad, mom, 1);
//     // auto exp_germline_mean = meiosis_diploid_mean_matrix(dad, mom);


//     // auto orig = f81::matrix(2*mu, freqs);
//     // auto exp_somatic_full = mitosis_diploid_matrix(orig);
//     // auto exp_somatic_nomut = mitosis_diploid_matrix(orig, 0);
//     // auto exp_somatic_posmut = exp_somatic_full - exp_somatic_nomut;
//     // auto exp_somatic_onemut = mitosis_diploid_matrix(orig, 1);
//     // auto exp_somatic_mean = mitosis_diploid_mean_matrix(orig);

//     // std::vector<TransitionMatrix> exp_full {{},{},exp_germline_full, exp_somatic_full, exp_somatic_full};
//     // std::vector<TransitionMatrix> exp_nomut {{},{},exp_germline_nomut, exp_somatic_nomut, exp_somatic_nomut};
//     // std::vector<TransitionMatrix> exp_posmut {{},{},exp_germline_posmut, exp_somatic_posmut, exp_somatic_posmut};
//     // std::vector<TransitionMatrix> exp_onemut {{},{},exp_germline_onemut, exp_somatic_onemut, exp_somatic_onemut};
//     // std::vector<TransitionMatrix> exp_mean {{},{},exp_germline_mean, exp_somatic_mean, exp_somatic_mean};

//     // FindMutations find_mutation {min_prob, r_graph, test_param_1};
//     // auto full_matrices = find_mutation.full_transition_matrices_;
//     // auto nomut_matrices = find_mutation.nomut_transition_matrices_;
//     // auto posmut_matrices = find_mutation.posmut_transition_matrices_;
//     // auto onemut_matrices = find_mutation.onemut_transition_matrices_;
//     // auto mean_matrices = find_mutation.mean_matrices_;
//     // for (int i = 0; i < 5; ++i) {
//     //     boost_check_matrix(exp_full[i], full_matrices[i]);
//     //     boost_check_matrix(exp_nomut[i], nomut_matrices[i]);
//     //     boost_check_matrix(exp_posmut[i], posmut_matrices[i]);
//     //     boost_check_matrix(exp_onemut[i], onemut_matrices[i]);
//     //     boost_check_matrix(exp_mean[i], mean_matrices[i]);
//     // }
// }


// BOOST_FIXTURE_TEST_CASE(test_operator, TrioWorkspace) {


// //     std::vector<peel::family_members_t> family {
// //             {1,4},
// //             {0,1,2},
// //             {0,3}
// //     };
// //     FindMutations::stats_t mutation_stats;
// // //    MutationStats stats(min_prob);
// //     FindMutations find_mutation{min_prob, r_graph, test_param_1};
// //     find_mutation(read_depths, ref_index, &mutation_stats);


// //     double scale = setup_workspace(ref_index, read_depths);

// //     //Test basic stats
// //     auto nomut_matrices = find_mutation.nomut_transition_matrices_;
// //     workspace.CleanupFast();
// //     dng::peel::up(workspace, family[0], nomut_matrices);
// //     dng::peel::to_father_fast(workspace, family[1], nomut_matrices);
// //     dng::peel::up(workspace, family[2], nomut_matrices);
// //     double result_nomut = log((workspace.lower[0] * workspace.upper[0]).sum());


// //     auto full_matrices = find_mutation.full_transition_matrices_;
// //     workspace.CleanupFast();
// //     dng::peel::up(workspace, family[0], full_matrices );
// //     dng::peel::to_father(workspace, family[1], full_matrices );
// //     dng::peel::up(workspace, family[2], full_matrices );
// //     double result_full = log((workspace.lower[0] * workspace.upper[0]).sum());


// //     // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
// //     const double expected_mup = -std::expm1(result_nomut - result_full);
// //     const double expected_lld = (result_full +scale) / M_LN10;
// //     const double expected_llh = result_full / M_LN10;

// //     BOOST_CHECK_CLOSE(expected_mup, mutation_stats.mup, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
// //     BOOST_CHECK_CLOSE(expected_lld, mutation_stats.lld, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
// // //    BOOST_CHECK_CLOSE(expected_llh, mutation_stats.llh, BOOST_CLOSE_PERCENTAGE_THRESHOLD);

// //     //Test posterior
// //     full_matrices = find_mutation.full_transition_matrices_;
// //     workspace.CleanupFast();
// //     dng::peel::up(workspace, family[0], full_matrices );
// //     dng::peel::to_father(workspace, family[1], full_matrices );
// //     dng::peel::up(workspace, family[2], full_matrices );

// //     dng::peel::up_reverse(workspace, family[2], full_matrices );
// //     dng::peel::to_father_reverse(workspace, family[1], full_matrices );
// //     dng::peel::up_reverse(workspace, family[0], full_matrices );

// //     dng::GenotypeArray expected_posterior;
// //     for (std::size_t i = 0; i < workspace.num_nodes; ++i) {
// //         expected_posterior = workspace.lower[i] * workspace.upper[i];
// //         expected_posterior /= expected_posterior.sum();

// //         boost_check_close_vector(expected_posterior,
// //         		mutation_stats.posterior_probabilities[i]);
// //     }



// }


// BOOST_FIXTURE_TEST_CASE(test_calculate_mutation_expected, TrioWorkspace) {

// //     FindMutations find_mutation{min_prob, r_graph, test_param_1};

// //     FindMutations::stats_t mutation_stats;
// // //    MutationStats mutation_stats(min_prob);
// //     find_mutation(read_depths, ref_index, &mutation_stats);

// //     BOOST_CHECK_SMALL(0.116189 - mutation_stats.mup, 1e-6);
// //     BOOST_CHECK_SMALL(-30.5967 - mutation_stats.lld, 1e-4);
// // //    BOOST_CHECK_SMALL(-6.68014 - mutation_stats.llh, 1e-5);
// //     BOOST_CHECK_SMALL(0.116189 - mutation_stats.mux, 1e-6);
// //     BOOST_CHECK_SMALL(0.116189 - mutation_stats.mu1p, 1e-6);

// //     BOOST_CHECK_EQUAL("GGxGG>GT", mutation_stats.dnt);
// //     BOOST_CHECK_EQUAL("LB-NA12878:Solexa-135852", mutation_stats.dnl);
// //     BOOST_CHECK_EQUAL(42, mutation_stats.dnq);

// // #if CALCULATE_ENTROPY == 1
// //     BOOST_CHECK_EQUAL(100, mutation_stats.dnc);
// // #endif

// }
