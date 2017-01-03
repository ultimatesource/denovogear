/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2016 Reed A. Cartwright
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
#ifndef DNG_FIND_MUTATIONS_ABSTRACT_H_
#define DNG_FIND_MUTATIONS_ABSTRACT_H_

#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <dng/relationship_graph.h>
#include <dng/genotyper.h>
#include <dng/mutation_stats.h>


namespace dng {

constexpr int SIZE4 = 4;


//TODO(SW): Eventually this will be just FindMutations.
//TODO(SW): The original FindMutations will become FindMutationsAutosomal
//TODO(SW): Once we refactor out task/loglike class, which can be a base class to this
class FindMutationsAbstract {


public:
    struct params_t {
        double theta;
        std::array<double, 4> nuc_freq;
        double ref_weight;

        dng::genotype::DirichletMultinomialMixture::params_t params_a;
        dng::genotype::DirichletMultinomialMixture::params_t params_b;
    };

    struct stats_t {
        float mup;
        float lld;
        [[deprecated]]float llh;
        float mux;

        bool has_single_mut;
        float mu1p;
        std::string dnt;
        std::string dnl;
        int32_t dnq;
        int32_t dnc;

        GenotypeArrayVector posterior_probabilities;
        GenotypeArrayVector genotype_likelihoods;
        std::vector<float> node_mup;
        std::vector<float> node_mu1p;
    };

    FindMutationsAbstract(double min_prob, const RelationshipGraph &graph,
            params_t params);

    virtual ~FindMutationsAbstract();

    virtual bool operator()(const pileup::RawDepths &depths, int ref_index,
            stats_t *stats) = 0;



protected:
    virtual void SetupTransitionMatrix() = 0;

    void SetupPopulationPriorDiploid();
    void SetupPopulationPriorHaploid();

    void Resize10To4(TransitionMatrix &matrix);
    void Resize10To4(GenotypeArray &array);

    bool CalculateMutationProb(MutationStats &mutation_stats);

    void CalculateDenovoMutation(MutationStats &mutation_stats);


    const dng::RelationshipGraph relationship_graph_;
    double min_prob_;
    params_t params_;
    dng::genotype::DirichletMultinomialMixture genotype_likelihood_;
    dng::peel::workspace_t work_full_;
    dng::peel::workspace_t work_nomut_;

    dng::GenotypeArray genotype_prior_[5]; // Holds P(G | theta)
    dng::TransitionMatrixVector full_transition_matrices_;
    dng::TransitionMatrixVector nomut_transition_matrices_;
    dng::TransitionMatrixVector posmut_transition_matrices_;
    dng::TransitionMatrixVector onemut_transition_matrices_;
    dng::TransitionMatrixVector mean_matrices_;

    std::vector<int> keep_library_index_;

#if CALCULATE_ENTROPY == 1
    std::array<double, 5> max_entropies_;
    std::vector<double> event_;
#endif
};
} // namespace dng

#endif /* DNG_FIND_MUTATIONS_ABSTRACT_H_ */
