/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2016,2017 Reed A. Cartwright
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

#include <boost/range/algorithm/transform.hpp>
#include <functional>
#include <iterator>

using namespace dng;

// helper function to simplify the code below
// Not optimized by any means, but reduces code mistakes
template<typename V>
V container_subtract(const V& a, const V& b) {
    V output;
    boost::transform(a, b, std::back_inserter(output),
        std::minus<typename V::value_type>());
    return output;
}

CallMutations::CallMutations(double min_prob, const RelationshipGraph &graph, params_t params)
        : LogProbability(graph, params), min_prob_{min_prob} {

    // Create Special Transition Matrices
    zero_mutation_matrices_ = CreateMutationMatrices(0);
    one_mutation_matrices_ = CreateMutationMatrices(1);
    mean_mutation_matrices_ = CreateMutationMatrices(MUTATIONS_MEAN);

    for(size_t j=0; j < oneplus_mutation_matrices_.size(); ++j) {
        oneplus_mutation_matrices_[j] = container_subtract(transition_matrices_[j],
            zero_mutation_matrices_[j]);
    }
}

// Returns true if a mutation was found and the record was modified
bool CallMutations::operator()(
    const pileup::allele_depths_t &depths, int num_alts, bool has_ref, stats_t *stats)
{
    // calculate genotype likelihoods and store in the lower library vector
    double scale = work_.SetGenotypeLikelihoods(genotyper_, depths, num_alts);

    // Set the prior probability of the founders given the reference
    work_.SetGermline(DiploidPrior(num_alts, has_ref), HaploidPrior(num_alts, has_ref));

    bool found = Calculate(stats, num_alts, has_ref);
    if(found && stats != nullptr) {
        stats->lld += scale/M_LN10;
    }

    return found;
}

bool CallMutations::Calculate(stats_t *stats, int num_alts, bool has_ref) {
    // Now peel numerator
    double numerator = graph_.PeelForwards(work_, zero_mutation_matrices_[num_alts]);

    // Calculate log P(Data ; model)
    double denominator = graph_.PeelForwards(work_, transition_matrices_[num_alts]);

    // Mutation Probability
    double mup = -std::expm1(numerator - denominator);

    if (mup < min_prob_) {
        return false;
    } else if(stats == nullptr) {
        return true;
    }
    stats->mup = mup;
    stats->lld = denominator/M_LN10;

    graph_.PeelBackwards(work_, transition_matrices_[num_alts]);

    // Genotype Likelihoods for Libraries
    size_t num_libraries = work_.library_nodes.second-work_.library_nodes.first;
    stats->genotype_likelihoods.resize(num_libraries);
    for (size_t u = 0; u < num_libraries; ++u) {
        size_t pos = work_.library_nodes.first + u;
        stats->genotype_likelihoods[u] = work_.lower[pos].log() / M_LN10;
    }

    // Posterior probabilities and best genotypes for all nodes
    stats->posterior_probabilities.resize(work_.num_nodes);
    stats->best_genotypes.resize(work_.num_nodes);
    stats->genotype_qualities.resize(work_.num_nodes);

    for (size_t i = 0; i < work_.num_nodes; ++i) {
        stats->posterior_probabilities[i] = (work_.upper[i] * work_.lower[i]);
        stats->posterior_probabilities[i] /= stats->posterior_probabilities[i].sum();

        size_t pos;
        double d = stats->posterior_probabilities[i].maxCoeff(&pos);
        stats->best_genotypes[i] = pos;
        stats->genotype_qualities[i] = dng::utility::lphred<int>(1.0 - d, 255);
    }

    // Expected Number of Mutations
    stats->mux = 0.0;
    for(size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        stats->mux += (work_.super[i] * (mean_mutation_matrices_[num_alts][i] *
                                      work_.lower[i].matrix()).array()).sum();
    }

    // Probability of at least 1 mutation at a node, given that there is at least 1 mutation in the graph
    stats->node_mup.resize(work_.num_nodes);
    for(size_t i = work_.founder_nodes.first; i < work_.founder_nodes.second; ++i) {
        stats->node_mup[i] = 0.0;
    }
    for (size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        stats->node_mup[i] = (work_.super[i] * (oneplus_mutation_matrices_[num_alts][i] *
                                          work_.lower[i].matrix()).array()).sum();
        stats->node_mup[i] /= mup;
    }

    // Probability of Exactly One Mutation
    // We will peel again, but this time with the zero matrix
    graph_.PeelForwards(work_, zero_mutation_matrices_[num_alts]);
    graph_.PeelBackwards(work_, zero_mutation_matrices_[num_alts]);

    double total = 0.0, max_coeff = -1.0;
    size_t dn_row = 0, dn_col = 0, dn_location = 0;

    stats->node_mu1p.resize(work_.num_nodes);
    for(size_t i = work_.founder_nodes.first; i < work_.founder_nodes.second; ++i) {
        stats->node_mu1p[i] = 0.0;
    }
    for (size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        work_.temp_buffer.resize(1,1);

        work_.temp_buffer = (work_.super[i].matrix() *
                                work_.lower[i].matrix().transpose()).array() *
                                one_mutation_matrices_[num_alts][i].array();
        std::size_t row, col;
        double temp = work_.temp_buffer.maxCoeff(&row, &col);
        if (temp > max_coeff) {
            max_coeff = temp;
            dn_row = row;
            dn_col = col;
            dn_location = i;
        }
        temp = work_.temp_buffer.sum();
        stats->node_mu1p[i] = temp;
        total += temp;
    }
    for (size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        stats->node_mu1p[i] /= total;
    }
    // total = P(1 mutation | D)/P(0 mutations | D) due to backwards algorithm
    stats->mu1p = total*exp(numerator-denominator);

    stats->dnq = dng::utility::lphred<int>(1.0 - (max_coeff / total), 255);
    stats->dnl = dn_location;
    stats->dnt_row = dn_row;
    stats->dnt_col = dn_col;
    return true;
}
