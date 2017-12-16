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

#include <numeric>

#include <dng/probability.h>
#include <dng/mutation.h>


using namespace dng;

// graph_ will be initialized before work_, so we can reference it.
LogProbability::LogProbability(RelationshipGraph graph, params_t params) :
    graph_{std::move(graph)},
    params_(std::move(params)),
    work_{graph_.CreateWorkspace()},
    genotyper_{params_.over_dispersion_hom, params_.over_dispersion_het, params_.ref_bias,
        params_.error_rate, params_.error_entropy}
{
    using namespace dng;

    // Create cache of population priors
    for(int i=0;i<diploid_prior_.size();++i) {
        diploid_prior_[i] = population_prior_diploid_ia(params_.theta,
            params_.asc_bias_hom, params_.asc_bias_het, i, true);
    }
    for(int i=0;i<diploid_prior_noref_.size();++i) {
        diploid_prior_noref_[i] = population_prior_diploid_ia(params_.theta,
            params_.asc_bias_hom, params_.asc_bias_het, i+1, false);
    }

    for(int i=0;i<haploid_prior_.size();++i) {
        haploid_prior_[i] = population_prior_haploid_ia(params_.theta,
            params_.asc_bias_hap, i, true);
    }
    for(int i=0;i<haploid_prior_noref_.size();++i) {
        haploid_prior_noref_[i] = population_prior_haploid_ia(params_.theta,
            params_.asc_bias_hap, i+1, false);
    }

    // Calculate mutation matrices
    transition_matrices_ = CreateMutationMatrices(MUTATIONS_ALL);

    // // Precalculate monomorphic histories (first 4 colors)
    // size_t num_libraries = work_.library_nodes.second - work_.library_nodes.first;
    // pileup::AlleleDepths depths{0,0,num_libraries,pileup::AlleleDepths::data_t(num_libraries, 0)};
    // for(int color=0; color<4;++color) {
    //     // setup monomorphic prior
    //     depths.color(color);
    //     GenotypeArray diploid_prior(1), haploid_prior(1);
    //     diploid_prior(0) = diploid_prior_[color](pileup::AlleleDepths::type_info_gt_table[color].indexes[0]);
    //     haploid_prior(0) = haploid_prior_[color](pileup::AlleleDepths::type_info_table[color].indexes[0]);
    //     work_.SetGermline(diploid_prior, haploid_prior);
    //     work_.SetGenotypeLikelihoods(genotyper_, depths);
    //     double logdata = graph_.PeelForwards(work_, transition_matrices_[color]);
    //     prob_monomorphic_[color] = exp(logdata);
    // }
}

// returns 'log10 P(Data ; model)-log10 scale' and log10 scaling.
LogProbability::value_t LogProbability::operator()(
    const pileup::allele_depths_t &depths, int num_alts, bool has_ref)
{
    // calculate genotype likelihoods and store in the lower library vector
    double scale = work_.SetGenotypeLikelihoods(genotyper_, depths, num_alts);

    // Set the prior probability of the founders given the reference
    work_.SetGermline(DiploidPrior(num_alts, has_ref), HaploidPrior(num_alts, has_ref));

    // Calculate log P(Data ; model)
    double logdata = graph_.PeelForwards(work_, transition_matrices_[num_alts]);

    return {logdata/M_LN10, scale/M_LN10};
}

// // Calculate the probability of a depths object considering only indexes
// // TODO: make indexes a property of the pedigree
// LogProbability::value_t LogProbability::operator()(const pileup::AlleleDepths &depths) {
//     const int ref_index = depths.type_info().reference;
//     const int color = depths.color();

//     double scale, logdata;
//     // For monomorphic sites we have pre-calculated the peeling part
//     if(color < 4) {
//         assert(depths.type_info().width == 1 && depths.type_gt_info().width == 1 && ref_index == color);
//         // Calculate genotype likelihoods
//         scale = work_.SetGenotypeLikelihoods(genotyper_, depths);

//         // Multiply our pre-calculated peeling results with the genotype likelihoods
//         logdata = prob_monomorphic_[color];
//         for(auto it = work_.lower.begin()+work_.library_nodes.first;
//             it != work_.lower.begin()+work_.library_nodes.second; ++it) {
//             logdata *= (*it)(0);
//         }
//         // convert to a log-likelihood
//         logdata = log(logdata);
//     } else {      
//         // Set the prior probability of the founders given the reference
//         GenotypeArray diploid_prior = DiploidPrior(ref_index, color);
//         GenotypeArray haploid_prior = HaploidPrior(ref_index, color);
//         work_.SetGermline(diploid_prior, haploid_prior);
  
//         // Calculate genotype likelihoods
//         scale = work_.SetGenotypeLikelihoods(genotyper_, depths);

//          // Calculate log P(Data ; model)
//         logdata = graph_.PeelForwards(work_, transition_matrices_[color]);
//     }
//     return {logdata/M_LN10, scale/M_LN10};
// }

// Construct the mutation matrices for each transition
TransitionMatrixVector dng::create_mutation_matrices(const RelationshipGraph &graph,
    int num_alleles, double entropy, const int mutype) {
    TransitionMatrixVector matrices(graph.num_nodes());
 
    for(size_t child = 0; child < graph.num_nodes(); ++child) {
        auto trans = graph.transition(child);
        if(trans.type == RelationshipGraph::TransitionType::Trio) {
            assert(graph.ploidy(child) == 2);
            auto dad = Mk::matrix(num_alleles, trans.length1, entropy);
            auto mom = Mk::matrix(num_alleles, trans.length2, entropy);
            matrices[child] = meiosis_matrix(graph.ploidy(trans.parent1), dad, graph.ploidy(trans.parent2), mom, mutype);
        } else if(trans.type == RelationshipGraph::TransitionType::Pair) {
            auto orig = Mk::matrix(num_alleles, trans.length1, entropy);
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
