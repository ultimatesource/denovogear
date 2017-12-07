/*
 * Copyright (c) 2016 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
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
#ifndef DNG_PROBABILITY_H
#define DNG_PROBABILITY_H

#include <array>
#include <vector>

#include <dng/genotyper.h>
#include <dng/relationship_graph.h>
#include <dng/mutation.h>

#include <dng/detail/unit_test.h>

namespace dng {

class LogProbability {
public:
    struct params_t {
        double theta;
        double hom_bias;
        double het_bias;
        double hap_bias;

        Genotyper::params_t genotyper_params;
    };

    struct value_t {
        double log_data;
        double log_scale;        
    };

    LogProbability(RelationshipGraph graph, params_t params);

    value_t operator()(const pileup::RawDepths &depths, int ref_index);
    value_t operator()(const pileup::AlleleDepths &depths);

    const peel::workspace_t& work() const { return work_; };

protected:
    using matrices_t = std::array<TransitionMatrixVector, pileup::AlleleDepths::type_info_table_length>;
    static constexpr int COLOR_ACGT = 40; // When color == 40 all four nucleotides are observed in order A,C,G,T

    matrices_t CreateMutationMatrices(const int mutype = MUTATIONS_ALL) const;

    GenotypeArray DiploidPrior(int num_alts, bool has_ref=true);
    GenotypeArray HaploidPrior(int num_alts, bool has_ref=true);

    RelationshipGraph graph_;
    params_t params_;
    peel::workspace_t work_; // must be declared after graph_ (see constructor)


    matrices_t transition_matrices_;

    double prob_monomorphic_[4];

    Genotyper genotyper_;

    GenotypeArray diploid_prior_[4]; // Holds P(G | theta)
    GenotypeArray haploid_prior_[4]; // Holds P(G | theta)
    GenotypeArray diploid_prior_noref_[4]; // Holds P(G | theta)
    GenotypeArray haploid_prior_noref_[4]; // Holds P(G | theta)

    DNG_UNIT_TEST_CLASS(unittest_dng_log_probability);
};

TransitionMatrixVector create_mutation_matrices(const RelationshipGraph &pedigree,
        int num_alleles, const int mutype = MUTATIONS_ALL);

TransitionMatrixVector create_mutation_matrices_subset(const TransitionMatrixVector &full_matrices, size_t color);

inline
LogProbability::matrices_t LogProbability::CreateMutationMatrices(const int mutype) const {
    matrices_t ret;
    // Construct the complete matrices
    auto full_matrix = create_mutation_matrices(graph_, params_.nuc_freq, mutype);

    // Extract relevant subsets of matrices
    for(size_t color = 0; color < dng::pileup::AlleleDepths::type_info_table_length; ++color) {
        ret[color] = create_mutation_matrices_subset(full_matrix, color);
    }
    return ret;
}

inline
GenotypeArray LogProbability::DiploidPrior(int num_alts, bool has_ref) {
    assert(num_alts >= 0);
    if(has_ref) {
        return (num_alts <= 3) ? diploid_prior_[num_alts]
            : population_prior_diploid_ia(params_.theta, params_.hom_bias, params_.het_bias, true);
    } else {
        assert(num_alts >= 1);
        return (num_alts < 5) ? diploid_prior_noref_[num_alts-1]
            : population_prior_diploid_ia(params_.theta, params_.hom_bias, params_.het_bias, false);
    }
}

inline
GenotypeArray LogProbability::HaploidPrior(int num_alts, bool has_ref) {
    assert(num_alts >= 0);
    if(has_ref) {
        return (num_alts <= 3) ? haploid_prior_[num_alts]
            : population_prior_haploid_ia(params_.theta, params_.hom_bias, params_.het_bias, true);
    } else {
        assert(num_alts >= 1);
        return (num_alts < 5) ? haploid_prior_noref_[num_alts-1]
            : population_prior_haploid_ia(params_.theta, params_.hom_bias, params_.het_bias, false);
    }
}

}; // namespace dng

#endif // DNG_PROBABILITY_H
