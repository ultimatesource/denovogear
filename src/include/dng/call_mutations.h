/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2016-2017 Reed A. Cartwright
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

class CallMutations : public LogProbability {
public:

    struct stats_t {
        double mup;
        double lld;
        double mux;

        double mu1p;
        std::string dnt;
        std::string dnl;
        int dnq;

        GenotypeArrayVector genotype_likelihoods;

        GenotypeArrayVector posterior_probabilities;
        std::vector<int> best_genotypes;
        std::vector<int> genotype_qualities;

        std::vector<double> node_mup;
        std::vector<double> node_mu1p;
    };

    CallMutations(double min_prob, const RelationshipGraph &graph,
            params_t params);

    bool operator()(const pileup::RawDepths &depths, int ref_index,
                    stats_t *stats);

    bool operator()(const pileup::AlleleDepths &depths,
                    stats_t *stats);
protected:
    bool Calculate(stats_t *stats, int color=COLOR_FULL);

    double min_prob_;

    matrices_t zero_mutation_matrices_;
    matrices_t one_mutation_matrices_;
    matrices_t oneplus_mutation_matrices_;
    matrices_t mean_mutation_matrices_;

    DNG_UNIT_TEST(test_constructor);
    DNG_UNIT_TEST(test_prior);
    DNG_UNIT_TEST(test_full_transition);
    DNG_UNIT_TEST(test_operator);

};
} // namespace dng

#endif /* DNG_FIND_MUTATIONS_H */
