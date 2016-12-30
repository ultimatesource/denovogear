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


#include <dng/call_mutations.h>

using namespace dng;

CallMutations::~CallMutations() {
	// TODO Auto-generated destructor stub
}

CallMutations::CallMutations(double min_prob, const RelationshipGraph &graph, params_t params)
        : LogProbability(graph, params), min_prob_{min_prob} {
    
    // Create Special Transition Matrices
    zero_mutation_matrices_ = CreateMutationMatrices(0);
    one_mutation_matrices_ = CreateMutationMatrices(1);
    mean_mutation_matrices_ = CreateMutationMatrices(MUTATIONS_MEAN);

    oneplus_mutation_matrices_.full.resize(transition_matrices_.full.size());
    for(int i=0; i < oneplus_mutation_matrices_.full.size(); ++i) {
        oneplus_mutation_matrices_.full[i] = transition_matrices_.full[i] - zero_mutation_matrices_.full[i];
    }
    for(size_t j=0; j < transition_matrices_.subsets[j].size(); ++j) {
        oneplus_mutation_matrices_.subsets[j].resize(transition_matrices_.subsets[j].size());
        for(int i=0; i < oneplus_mutation_matrices_.subsets[j].size(); ++i) {
            oneplus_mutation_matrices_.subsets[j][i] = transition_matrices_.subsets[j][i] - zero_mutation_matrices_.subsets[j][i];
        }        
    }
}

// Returns true if a mutation was found and the record was modified
bool CallMutations::operator()(const RawDepths &depths,
        int ref_index, stats_t *stats) {
    assert(stats != nullptr);

    // Genotype Likelihoods
    double scale = work_.SetGenotypeLikelihoods(genotyper_, depths, ref_index);

    // Set the prior probability of the founders given the reference
    work_.SetFounders(diploid_prior_[ref_index], haploid_prior_[ref_index]);

    // Calculate log P(Data ; model)
    double denominator = pedigree_.PeelForwards(work_, transition_matrices_.full);

    // Now peel numerator
    double numerator = pedigree_.PeelForwards(work_, zero_mutation_matrices_.full);

    // Mutation Probability
    double mup = -std::expm1(numerator - denominator);

    if (mup < min_prob_) {
        return false;
    }
    stats->mup = mup;
    return true;
}
