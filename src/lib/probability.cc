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

    mutation::population_prior_check(params_.theta, params_.ref_bias_hom, params_.ref_bias_het,
        params_.ref_bias_hap, params_.k_alleles);

    // Create cache of population priors
    for(int i=0;i<diploid_prior_.size();++i) {
        diploid_prior_[i] = mutation::population_prior_diploid(i+1, params_.theta,
            params_.ref_bias_hom, params_.ref_bias_het, params_.k_alleles);
    }

    for(int i=0;i<haploid_prior_.size();++i) {
        haploid_prior_[i] = mutation::population_prior_haploid(i+1, params_.theta,
            params_.ref_bias_hap, params_.k_alleles);
    }

    // Calculate mutation matrices
    transition_matrices_ = CreateMutationMatrices(mutation::transition_t{});

    // Precalculate monomorphic histories
    size_t num_libraries = work_.library_nodes.second - work_.library_nodes.first;
    work_.matrix_index = 0;
    work_.ClearGenotypeLikelihoods(1);
    work_.SetGermline(DiploidPrior(1), HaploidPrior(1));
    ln_monomorphic_ = graph_.PeelForwards(work_, transition_matrices_[0]);
}

