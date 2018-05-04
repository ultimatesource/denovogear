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
    CallMutations(const RelationshipGraph &graph, params_t params);

    struct stats_t {
        double mutq;
        double mutx;
        double lld;
        //double lld1;

        double dnp;
        int dnt_row;
        int dnt_col;
        int dnl;
        int dnq;

        double quality;

        GenotypeArrayVector genotype_likelihoods;

        GenotypeArrayVector posterior_probabilities;
        std::vector<int> best_genotypes;
        std::vector<int> genotype_qualities;

        std::vector<double> node_mutp;
        std::vector<double> node_dnp;

        double ln_zero;
        double ln_all;
        bool has_single_mut;
    };

    double PeelNoMutations();

    bool CalculateMutationStats(genotype::Mode mode, stats_t *stats);

    Probability::logdiff_t CalculateMUQ(stats_t *stats);

    double CalculateSingleMutationStats(stats_t *stats);

    bool all_variants() const { return all_variants_; }
    double quality_threshold() const { return min_quality_; }

    void all_variants(bool all) { all_variants_ = all; }
    void quality_threshold(double min_quality) { min_quality_ = min_quality; }
    void quality_threshold(double min_quality, bool all ) {
        min_quality_ = min_quality;
        all_variants_ = all;
    }

protected:

    double min_quality_{0};
    bool all_variants_{false};

    matrices_t zero_mutation_matrices_;
    matrices_t one_mutation_matrices_;
    matrices_t oneplus_mutation_matrices_;
    matrices_t mean_mutation_matrices_;

    std::array<double, MAXIMUM_NUMBER_ALLELES> one_mutation_prior_; 

    DNG_UNIT_TEST_CLASS(unittest_dng_call_mutations);
};

inline
double CallMutations::PeelNoMutations() {
    const int matrix_index = work_.matrix_index;
    return graph_.PeelForwards(work_, zero_mutation_matrices_[matrix_index]);
}

} // namespace dng

#endif /* DNG_FIND_MUTATIONS_H */
