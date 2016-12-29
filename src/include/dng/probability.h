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

#include <dng/likelihood.h>
#include <dng/relationship_graph.h>
#include <dng/mutation.h>

namespace dng {

class LogProbability {
public:
    struct params_t {
        double theta;
        std::array<double, 4> nuc_freq;
        double ref_weight;

        dng::genotype::DirichletMultinomialMixture::params_t params_a;
        dng::genotype::DirichletMultinomialMixture::params_t params_b;
    };

    struct stats_t {
        double log_data;
        double log_scale;        
    };

    LogProbability(RelationshipGraph pedigree, params_t params);

    stats_t operator()(const std::vector<depth_t> &depths, int ref_index);
    stats_t operator()(const pileup::AlleleDepths &depths, const std::vector<size_t>& indexes);

protected:
    struct matrices_t {
        TransitionVector full;
        TransitionVector subsets[pileup::AlleleDepths::type_info_table_length];
    };

    matrices_t CreateMutationMatrices(const int mutype = MUTATIONS_ALL) const;

    double CalculateGenotypeLikelihoods(const pileup::AlleleDepths &depths, const std::vector<size_t>& indexes);

    RelationshipGraph pedigree_;
    params_t params_;
    peel::workspace_t work_; // must be declared after pedigree_ (see constructor)


    matrices_t transition_matrices_;

    double prob_monomorphic_[4];

    // Model genotype likelihoods as a mixture of two dirichlet multinomials
    // TODO: control these with parameters
    genotype::DirichletMultinomialMixture genotype_likelihood_;

    GenotypeArray diploid_prior_[5]; // Holds P(G | theta)
    GenotypeArray haploid_prior_[5]; // Holds P(G | theta)
};

TransitionVector create_mutation_matrices(const RelationshipGraph &pedigree,
        const std::array<double, 4> &nuc_freq, const int mutype = MUTATIONS_ALL);

TransitionVector create_mutation_matrices_subset(const TransitionVector &full_matrices, size_t color);

inline
LogProbability::matrices_t LogProbability::CreateMutationMatrices(const int mutype) const {
    matrices_t ret;
    // Construct the complete matrices
    ret.full = create_mutation_matrices(pedigree_, params_.nuc_freq, mutype);

    // Extract relevant subsets of matrices
    for(size_t color = 0; color < dng::pileup::AlleleDepths::type_info_table_length; ++color) {
        ret.subsets[color] = create_mutation_matrices_subset(ret.full, color);
    }
    return ret;
}

}; // namespace dng

#endif // DNG_PROBABILITY_H
