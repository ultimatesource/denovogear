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

// pedigree_ will be initialized before work_, so we will reference it.
LogProbability::LogProbability(RelationshipGraph pedigree, params_t params) :
    pedigree_{std::move(pedigree)},
    params_(std::move(params)),
    genotype_likelihood_{params.params_a, params.params_b},
    work_{pedigree_.CreateWorkspace()} {

    using namespace dng;

    // Use a parent-independent mutation model, which produces a
    // beta-binomial
    genotype_prior_[0] = population_prior(params_.theta, params_.nuc_freq, {params_.ref_weight, 0, 0, 0});
    genotype_prior_[1] = population_prior(params_.theta, params_.nuc_freq, {0, params_.ref_weight, 0, 0});
    genotype_prior_[2] = population_prior(params_.theta, params_.nuc_freq, {0, 0, params_.ref_weight, 0});
    genotype_prior_[3] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, params_.ref_weight});
    genotype_prior_[4] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, 0});

    // Calculate mutation matrices
    full_transition_matrices_.assign(work_.num_nodes, {});
 
    for(size_t child = 0; child < work_.num_nodes; ++child) {
        auto trans = pedigree_.transitions()[child];
        if(trans.type == RelationshipGraph::TransitionType::Germline) {
            auto dad = f81::matrix(trans.length1, params_.nuc_freq);
            auto mom = f81::matrix(trans.length2, params_.nuc_freq);

            full_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom);
        } else if(trans.type == RelationshipGraph::TransitionType::Somatic ||
                  trans.type == RelationshipGraph::TransitionType::Library) {
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);

            full_transition_matrices_[child] = mitosis_diploid_matrix(orig);
        } else {
            full_transition_matrices_[child] = {};
        }
    }

    // Extract relevant subsets of matrices
    for(size_t color = 0; color < dng::pileup::AlleleDepths::type_info_table_length; ++color) {
        int width = dng::pileup::AlleleDepths::type_info_gt_table[color].width;
        // Resize our subsets to the right width
        transition_matrices_[color].resize(full_transition_matrices_.size());
        // enumerate over all children
        for(size_t child = 0; child < work_.num_nodes; ++child) {
            auto trans = pedigree_.transitions()[child];
            if(trans.type == RelationshipGraph::TransitionType::Germline) {
                // resize transition matrix to w*w,w
                transition_matrices_[color][child].resize(width*width,width);
                // Assume column major order which is the default
                for(int a = 0; a < width; ++a) {
                    int ga = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[a];
                    for(int b = 0; b < width; ++b) {
                        int gb = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[b];
                        // copy correct value from the full matrix to the subset matrix
                        int x = a*width+b;
                        int gx = ga*10+gb;
                        for(int y = 0; y < width; ++y) {
                            int gy = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[y];
                            // copy correct value from the full matrix to the subset matrix
                            transition_matrices_[color][child](x,y) = full_transition_matrices_[child](gx,gy); 
                        }
                    }
                }
            } else if(trans.type == RelationshipGraph::TransitionType::Somatic ||
                      trans.type == RelationshipGraph::TransitionType::Library) {
                transition_matrices_[color][child].resize(width,width);
                // Assume column major order which is the default
                for(int x = 0; x < width; ++x) {
                    int gx = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[x];
                    for(int y = 0; y < width; ++y) {
                        int gy = dng::pileup::AlleleDepths::type_info_gt_table[color].indexes[y];
                        // copy correct value from the full matrix to the subset matrix
                        transition_matrices_[color][child](x,y) = full_transition_matrices_[child](gx,gy); 
                    }
                }
            } else {
                transition_matrices_[color][child] = {};
            }
        }
    }

    // Precalculate monomorphic histories (first 4 colors)
    size_t num_libraries = work_.library_nodes.second - work_.library_nodes.first;
    pileup::AlleleDepths depths{0,0,num_libraries,pileup::AlleleDepths::data_t(num_libraries, 0)};
    std::vector<size_t> indexes(num_libraries,0);
    std::iota(indexes.begin(),indexes.end(),0);

    for(int color=0; color<4;++color) {
        // setup monomorphic prior
        depths.color(color);
        GenotypeArray prior(1);
        prior(0) = genotype_prior_[color](pileup::AlleleDepths::type_info_gt_table[color].indexes[0]);
        work_.SetFounders(prior);
        double scale = genotype_likelihood_(depths, indexes, work_.lower.begin()+work_.library_nodes.first);
        double logdata = pedigree_.PeelForwards(work_, transition_matrices_[color]);
        prob_monomorphic_[color] = exp(logdata);
    }
}

// returns 'log10 P(Data ; model)-log10 scale' and log10 scaling.
LogProbability::stats_t LogProbability::operator()(const std::vector<depth_t> &depths,
                               int ref_index) {
    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::tie(work_.lower[work_.library_nodes.first + u], stemp) =
            genotype_likelihood_(depths[u], ref_index);
        scale += stemp;
    }

    // Set the prior probability of the founders given the reference
    work_.SetFounders(genotype_prior_[ref_index]);

    // Calculate log P(Data ; model)
    double logdata = pedigree_.PeelForwards(work_, full_transition_matrices_);

    return {logdata/M_LN10, scale/M_LN10};
}

// Calculate the probability of a depths object considering only indexes
// TODO: make indexes a property of the pedigree
LogProbability::stats_t LogProbability::operator()(const pileup::AlleleDepths &depths, const std::vector<size_t>& indexes) {
    int ref_index = depths.type_info().reference;
    int color = depths.color();
    size_t gt_width = depths.type_gt_info().width;

    double scale, logdata;
    // For monomorphic sites we have pre-calculated the peeling part
    if(color < 4) {
        assert(gt_width == 1 && ref_index == color);
        // Calculate genotype likelihoods
        scale = genotype_likelihood_(depths, indexes, work_.lower.begin()+work_.library_nodes.first);

        // Multiply our pre-calculated peeling results with the genotype likelihoods
        logdata = prob_monomorphic_[color];
        for(auto it = work_.lower.begin()+work_.library_nodes.first;
            it != work_.lower.begin()+work_.library_nodes.second; ++it) {
            logdata *= (*it)(0);
        }
        // convert to a log-likelihood
        logdata = log(logdata);
    } else {
        // Set the prior probability of the founders given the reference
        GenotypeArray prior(gt_width);
        for(int i=0;i<gt_width;++i) {
            prior(i) = genotype_prior_[ref_index](depths.type_gt_info().indexes[i]);
        }
        work_.SetFounders(prior);
    
        // Calculate genotype likelihoods
        scale = genotype_likelihood_(depths, indexes, work_.lower.begin()+work_.library_nodes.first);

         // Calculate log P(Data ; model)
        logdata = pedigree_.PeelForwards(work_, transition_matrices_[color]);
    }
    return {logdata/M_LN10, scale/M_LN10};
}
