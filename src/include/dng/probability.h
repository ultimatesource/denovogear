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
        double ref_bias_hom;
        double ref_bias_het;
        double ref_bias_hap;
        
        double over_dispersion_hom;
        double over_dispersion_het;
        double sequencing_bias;
        double error_rate;
        double lib_k_alleles;

        double k_alleles;
    };

    struct value_t {
        double log_data;
        double log_scale;        
    };

    LogProbability(RelationshipGraph graph, params_t params);

    template<typename A>
    value_t CalculateLLD(const A &depths, int num_alts, bool has_ref=true);

    template<typename A>
    value_t operator()(const A &depths, int num_alts, bool has_ref=true) {
        return CalculateLLD(depths, num_alts, has_ref);
    }

    const peel::workspace_t& work() const { return work_; };

protected:
    using matrices_t = std::array<TransitionMatrixVector, 4>;

    matrices_t CreateMutationMatrices(const int mutype = MUTATIONS_ALL) const;

    GenotypeArray DiploidPrior(int num_obs_alleles, bool has_ref);
    GenotypeArray HaploidPrior(int num_obs_alleles, bool has_ref);

    RelationshipGraph graph_;
    params_t params_;
    peel::workspace_t work_; // must be declared after graph_ (see constructor)

    matrices_t transition_matrices_;

    double prob_monomorphic_;

    Genotyper genotyper_;

    std::array<GenotypeArray,4> diploid_prior_; // Holds P(G | theta)
    std::array<GenotypeArray,4> haploid_prior_; // Holds P(G | theta)
    
    std::array<GenotypeArray,4> diploid_prior_noref_; // Holds P(G | theta)
    std::array<GenotypeArray,4> haploid_prior_noref_; // Holds P(G | theta)

    DNG_UNIT_TEST_CLASS(unittest_dng_log_probability);
};

// returns 'log10 P(Data ; model)-log10 scale' and log10 scaling.
template<typename A>
LogProbability::value_t LogProbability::CalculateLLD(
    const A &depths, int num_observed_alleles, bool has_ref)
{
    if(num_alts >= transition_matrices_.size()) {
        num_alts = transition_matrices_.size()-1;
    }

    // calculate genotype likelihoods and store in the lower library vector
    double scale = work_.SetGenotypeLikelihoods(genotyper_, depths, num_alts);
    double logdata;

    if(num_alts == 0) {
        // Use cached value for monomorphic sites instead of peeling.
        logdata = prob_monomorphic_;
        for(auto it = work_.lower.begin()+work_.library_nodes.first;
            it != work_.lower.begin()+work_.library_nodes.second; ++it) {
            logdata *= (*it)(0);
        }
        // convert to a log-likelihood
        logdata = log(logdata);
    } else {
        // Set the prior probability of the founders given the reference
        work_.SetGermline(DiploidPrior(num_alts, has_ref), HaploidPrior(num_alts, has_ref));

        // Calculate log P(Data ; model)
        logdata = graph_.PeelForwards(work_, transition_matrices_[num_alts]);
    }
    return {logdata/M_LN10, scale/M_LN10};
}

TransitionMatrixVector create_mutation_matrices(const RelationshipGraph &pedigree,
        int num_alleles, double num_mutants, const int mutype = MUTATIONS_ALL);

inline
LogProbability::matrices_t LogProbability::CreateMutationMatrices(const int mutype) const {
    // Construct the complete matrices
    matrices_t ret;
    for(int i=0;i<ret.size();++i) {
        ret[i] = create_mutation_matrices(graph_, i+1, params_.k_alleles, mutype);
    }
    return ret;
}

inline
GenotypeArray LogProbability::DiploidPrior(int num_obs_alleles, bool has_ref) {
    assert(num_obs_alleles >= 0);
    if(has_ref) {
        return (num_obs_alleles < diploid_prior_.size()) ? diploid_prior_[num_obs_alleles]
            : population_prior_diploid(num_obs_alleles, params_.theta, params_.ref_bias_hom,
                params_.ref_bias_het, params_.k_alleles, true);
    } else {
        assert(num_obs_alleles >= 1);
        return (num_obs_alleles-1 < diploid_prior_noref_.size()) ? diploid_prior_noref_[num_obs_alleles-1]
            : population_prior_diploid(num_obs_alleles, params_.theta, params_.ref_bias_hom,
                 params_.ref_bias_het, params_.k_alleles, false);
    }
}

inline
GenotypeArray LogProbability::HaploidPrior(int num_obs_alleles, bool has_ref) {
    assert(num_obs_alleles >= 0);
    if(has_ref) {
        return (num_obs_alleles < haploid_prior_.size()) ? haploid_prior_[num_obs_alleles]
            : population_prior_haploid(num_obs_alleles, params_.theta, params_.ref_bias_hap, params_.k_alleles, true);
    } else {
        assert(num_obs_alleles >= 1);
        return (num_obs_alleles-1 < haploid_prior_noref_.size()) ? haploid_prior_noref_[num_obs_alleles-1]
            : population_prior_haploid(num_obs_alleles, params_.theta, params_.ref_bias_hap, params_.k_alleles, false);
    }
}

template<typename A>
inline
LogProbability::params_t get_model_parameters(const A& a) {
    LogProbability::params_t ret;
    
    ret.theta = a.theta;
    ret.ref_bias_hom = a.ref_bias_hom;
    ret.ref_bias_het = a.ref_bias_het;
    ret.ref_bias_hap = a.ref_bias_hap;
        
    ret.over_dispersion_hom = a.lib_overdisp_hom;
    ret.over_dispersion_het = a.lib_overdisp_het;
    ret.sequencing_bias = a.lib_bias;
    ret.error_rate = a.lib_error;
    ret.lib_k_alleles = a.lib_kbases;

    ret.k_alleles = a.kalleles;

    return ret;
}

}; // namespace dng

#endif // DNG_PROBABILITY_H
