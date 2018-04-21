/*
 * Copyright (c) 2016-2018 Reed A. Cartwright
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

#include <numeric>

#include <dng/probability.h>
#include <dng/mutation.h>

using namespace dng;

// graph_ will be initialized before work_, so we can reference it.
Probability::Probability(RelationshipGraph graph, params_t params) :
    graph_{std::move(graph)},
    params_(std::move(params)),
    work_{graph_.CreateWorkspace()},
    genotyper_{params_.over_dispersion_hom, params_.over_dispersion_het, params_.sequencing_bias,
        params_.error_rate, params_.lib_k_alleles}
{
    using namespace dng;

    // Create cache of population priors
    for(int i=0;i<diploid_prior_.size();++i) {
        diploid_prior_[i] = population_prior_diploid(i+1, params_.theta,
            params_.ref_bias_hom, params_.ref_bias_het, params_.k_alleles, true);
    }
    for(int i=0;i<diploid_prior_noref_.size();++i) {
        diploid_prior_noref_[i] = population_prior_diploid(i+2, params_.theta,
            params_.ref_bias_hom, params_.ref_bias_het, params_.k_alleles, false);
    }

    for(int i=0;i<haploid_prior_.size();++i) {
        haploid_prior_[i] = population_prior_haploid(i+1, params_.theta,
            params_.ref_bias_hap, params_.k_alleles, true);
    }
    for(int i=0;i<haploid_prior_noref_.size();++i) {
        haploid_prior_noref_[i] = population_prior_haploid(i+2, params_.theta,
            params_.ref_bias_hap, params_.k_alleles, false);
    }

    // Calculate mutation matrices
    transition_matrices_ = CreateMutationMatrices(MUTATIONS_ALL);

    // Precalculate monomorphic histories
    size_t num_libraries = work_.library_nodes.second - work_.library_nodes.first;
    work_.SetGermline(diploid_prior_[0], haploid_prior_[0]);
    boost::multi_array<int32_t, 2> genotype_likelihoods(utility::make_array(num_libraries, 1u));
    for(auto &&a : genotype_likelihoods) {
        a[0] = 1.0;
    }
    work_.CopyGenotypeLikelihoods(genotype_likelihoods);
    double logdata = graph_.PeelForwards(work_, transition_matrices_[0]);
    prob_monomorphic_ = exp(logdata);
}

// Construct the mutation matrices for each transition
TransitionMatrixVector dng::create_mutation_matrices(const RelationshipGraph &graph,
    int num_obs_alleles, double k_alleles, const int mutype) {
    TransitionMatrixVector matrices(graph.num_nodes());
 
    for(size_t child = 0; child < graph.num_nodes(); ++child) {
        auto trans = graph.transition(child);
        if(trans.type == RelationshipGraph::TransitionType::Trio) {
            assert(graph.ploidy(child) == 2);
            auto dad = Mk::matrix(num_obs_alleles, trans.length1, k_alleles);
            auto mom = Mk::matrix(num_obs_alleles, trans.length2, k_alleles);
            matrices[child] = meiosis_matrix(graph.ploidy(trans.parent1),
                dad, graph.ploidy(trans.parent2), mom, mutype);
        } else if(trans.type == RelationshipGraph::TransitionType::Pair) {
            auto orig = Mk::matrix(num_obs_alleles, trans.length1, k_alleles);
            if(graph.ploidy(child) == 1) {
                matrices[child] = gamete_matrix(graph.ploidy(trans.parent1), orig, mutype);
            } else {
                assert(graph.ploidy(child) == 2);
                matrices[child] = mitosis_matrix(graph.ploidy(trans.parent1), orig, mutype);
            }
        } else {
            matrices[child] = {};
        }
    }
    return matrices;
}
