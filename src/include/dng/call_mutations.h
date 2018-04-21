/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2016-2018 Reed A. Cartwright
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

#pragma once
#ifndef DNG_CALL_MUTATIONS_H
#define DNG_CALL_MUTATIONS_H

#include <dng/probability.h>
#include <dng/depths.h>

#include <dng/detail/unit_test.h>

namespace dng {

class CallMutations : public Probability {
public:
    CallMutations(double min_prob, const RelationshipGraph &graph, params_t params);

    struct stats_t {
        double ln_zero;
        double ln_all;

        double mup;
        double lld;
        double lld1;
        double mux;

        double mu1p;
        int dnt_row;
        int dnt_col;
        int dnl;
        int dnq;

        double quality;

        GenotypeArrayVector genotype_likelihoods;

        GenotypeArrayVector posterior_probabilities;
        std::vector<int> best_genotypes;
        std::vector<int> genotype_qualities;

        std::vector<double> node_mup;
        std::vector<double> node_mu1p;
    };

    bool CalculateMutationStats(stats_t *stats);

    double CalculateMUP(stats_t *stats);

    double CalculateMU1P(stats_t *stats);

protected:

    double min_prob_;

    matrices_t zero_mutation_matrices_;
    matrices_t one_mutation_matrices_;
    matrices_t oneplus_mutation_matrices_;
    matrices_t mean_mutation_matrices_;

    std::array<double, MAXIMUM_NUMBER_ALLELES> log10_one_mutation_; 

    DNG_UNIT_TEST_CLASS(unittest_dng_call_mutations);
};

inline
double CallMutations::CalculateMUP(stats_t *stats) {
    const int matrix_index = work_.matrix_index;
    // Now peel numerator
    double numerator = graph_.PeelForwards(work_, zero_mutation_matrices_[matrix_index]);
    // Calculate log P(Data ; model)
    double denominator = graph_.PeelForwards(work_, transition_matrices_[matrix_index]);
    // Mutation Probability
    double mup = -std::expm1(numerator - denominator);
    if(stats != nullptr) {
        stats->mup = mup;
        stats->ln_all = denominator;
        stats->ln_zero = numerator;
    }
    return mup;
}


} // namespace dng

#endif /* DNG_FIND_MUTATIONS_H */
