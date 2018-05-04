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

class Probability {
public:
    struct params_t;
    struct logdiff_t;
    static constexpr int MAXIMUM_NUMBER_ALLELES{4};

    static int adjust_num_obs_alleles(int num) {
        assert(num >= 1);
        return (num <= MAXIMUM_NUMBER_ALLELES) ? num : MAXIMUM_NUMBER_ALLELES;
    }

    Probability(RelationshipGraph graph, params_t params);

    template<typename A>
    void SetupWorkspace(const A &depths, int num_obs_alleles, genotype::Mode mode);

    double CalculateLLD();

    double PeelOnlyReference(genotype::Mode mode);

    logdiff_t CalculateMONO(genotype::Mode mode);

    template<typename A>
    double CalculateLLD(const A &depths, int num_obs_alleles);

    const peel::workspace_t& work() const { return work_; }

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

protected:
    using matrices_t = std::array<TransitionMatrixVector, MAXIMUM_NUMBER_ALLELES>;

    matrices_t CreateMutationMatrices(const int mutype = MUTATIONS_ALL) const;

    GenotypeArray DiploidPrior(int num_obs_alleles, bool has_ref);
    GenotypeArray HaploidPrior(int num_obs_alleles, bool has_ref);

    RelationshipGraph graph_;
    params_t params_;
    peel::workspace_t work_; // must be declared after graph_ (see constructor)

    matrices_t transition_matrices_;

    double ln_monomorphic_;

    Genotyper genotyper_;

    using prior_t = std::array<GenotypeArray, MAXIMUM_NUMBER_ALLELES>;

    prior_t diploid_prior_; // Holds P(G | theta)
    prior_t haploid_prior_; // Holds P(G | theta)
    
    prior_t diploid_prior_noref_; // Holds P(G | theta)
    prior_t haploid_prior_noref_; // Holds P(G | theta)

    DNG_UNIT_TEST_CLASS(unittest_dng_log_probability);
};

struct Probability::logdiff_t {
    double left;
    double right;

    double value() const {
        assert(left <= right);
        return left - right;
    }

    double phred_score() const { return (-10.0/M_LN10)*value(); }
    double prob() const { return exp(value()); }
    double not_prob() const { return -expm1(value()); }
};

template<typename A>
void Probability::SetupWorkspace(const A &depths, int num_obs_alleles, genotype::Mode mode) {
    assert(num_obs_alleles >= 1);
    num_obs_alleles = adjust_num_obs_alleles(num_obs_alleles);
    work_.matrix_index = num_obs_alleles-1;

    work_.CalculateGenotypeLikelihoods(genotyper_, depths, num_obs_alleles, mode);
    work_.SetGermline(DiploidPrior(num_obs_alleles, true), HaploidPrior(num_obs_alleles, true));
}

inline
Probability::logdiff_t Probability::CalculateMONO(genotype::Mode mode) {
    assert(work_.matrix_index >= 0 && work_.matrix_index < transition_matrices_.size());
    if(work_.matrix_index == 0) {
        if(mode != genotype::Mode::Likelihood) {
            work_.ExpGenotypeLikelihoods();
        }
        return {ln_monomorphic_, ln_monomorphic_};
    }
    double ln_mono = PeelOnlyReference(mode);
    if(mode != genotype::Mode::Likelihood) {
        work_.ExpGenotypeLikelihoods();
    }
    double ln_data = graph_.PeelForwards(work_, transition_matrices_[work_.matrix_index]);
    return {ln_mono, ln_data};
}

inline
double Probability::PeelOnlyReference(genotype::Mode mode) {
    // Use cached value for monomorphic sites instead of peeling.
    if(work_.matrix_index == 0) {
        return ln_monomorphic_;
    }
    double ln_mono = ln_monomorphic_;
    for(auto &&lower : work_.make_library_slice(work_.lower)) {
        if(mode == genotype::Mode::Likelihood) {
            ln_mono += log(lower(0));
        } else {
            ln_mono += lower(0);
        }
    }
    return ln_mono;
}

// returns 'log10 P(Data ; model)'
inline
double Probability::CalculateLLD() {
    double ln_data = graph_.PeelForwards(work_, transition_matrices_[work_.matrix_index]);
    return (ln_data+work_.ln_scale)/M_LN10;
}

// returns 'log10 P(Data ; model)'
template<typename A>
double Probability::CalculateLLD(const A &depths, int num_obs_alleles)
{
    num_obs_alleles = adjust_num_obs_alleles(num_obs_alleles);
    work_.matrix_index = num_obs_alleles-1;

    // calculate genotype likelihoods and store in the lower library vector
    if(num_obs_alleles > 1) {
        SetupWorkspace(depths, num_obs_alleles, genotype::Mode::Likelihood);
        return CalculateLLD();
    }
    // Use cached value for monomorphic sites instead of peeling.
    work_.CalculateGenotypeLikelihoods(genotyper_, depths, 1, genotype::Mode::Likelihood);
    return (ln_monomorphic_ + work_.ln_scale)/M_LN10;
}

TransitionMatrixVector create_mutation_matrices(const RelationshipGraph &pedigree,
        int num_alleles, double num_mutants, const int mutype = MUTATIONS_ALL);

inline
Probability::matrices_t Probability::CreateMutationMatrices(const int mutype) const {
    // Construct the complete matrices
    matrices_t ret;
    for(int i=0;i<ret.size();++i) {
        ret[i] = create_mutation_matrices(graph_, i+1, params_.k_alleles, mutype);
    }
    return ret;
}

inline
GenotypeArray Probability::DiploidPrior(int num_obs_alleles, bool has_ref) {
    assert(num_obs_alleles >= 1);
    if(has_ref) {
        return (num_obs_alleles-1 < diploid_prior_.size()) ? diploid_prior_[num_obs_alleles-1]
            : population_prior_diploid(num_obs_alleles, params_.theta, params_.ref_bias_hom,
                params_.ref_bias_het, params_.k_alleles, true);
    } else {
        assert(num_obs_alleles >= 2);
        return (num_obs_alleles-2 < diploid_prior_noref_.size()) ? diploid_prior_noref_[num_obs_alleles-2]
            : population_prior_diploid(num_obs_alleles, params_.theta, params_.ref_bias_hom,
                 params_.ref_bias_het, params_.k_alleles, false);
    }
}

inline
GenotypeArray Probability::HaploidPrior(int num_obs_alleles, bool has_ref) {
    assert(num_obs_alleles >= 1);
    if(has_ref) {
        return (num_obs_alleles < haploid_prior_.size()) ? haploid_prior_[num_obs_alleles-1]
            : population_prior_haploid(num_obs_alleles, params_.theta, params_.ref_bias_hap, params_.k_alleles, true);
    } else {
        assert(num_obs_alleles >= 2);
        return (num_obs_alleles-1 < haploid_prior_noref_.size()) ? haploid_prior_noref_[num_obs_alleles-2]
            : population_prior_haploid(num_obs_alleles, params_.theta, params_.ref_bias_hap, params_.k_alleles, false);
    }
}

template<typename A>
inline
Probability::params_t get_model_parameters(const A& a) {
    Probability::params_t ret;
    
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
