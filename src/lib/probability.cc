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
    genotyper_{params.genotyper_params}
    {

    using namespace dng;

    // Use a parent-independent mutation model, which produces a beta-binomial
    for(int i=0; i < 5; i++) {
        std::array<double, 4> weights = {{0,0,0,0}};
        if(i < 4) {
            weights[i] = params_.ref_weight;
        }
        haploid_prior_[i] = population_prior_haploid(params_.theta, params_.nuc_freq, weights);
        diploid_prior_[i] = population_prior_diploid(params_.theta, params_.nuc_freq, weights);
    }

    // Calculate mutation matrices
    transition_matrices_ = CreateMutationMatrices(MUTATIONS_ALL);

    // Precalculate monomorphic histories (first 4 colors)
    size_t num_libraries = work_.library_nodes.second - work_.library_nodes.first;
    pileup::AlleleDepths depths{0,0,num_libraries,pileup::AlleleDepths::data_t(num_libraries, 0)};
    for(int color=0; color<4;++color) {
        // setup monomorphic prior
        depths.color(color);
        GenotypeArray diploid_prior(1), haploid_prior(1);
        diploid_prior(0) = diploid_prior_[color](pileup::AlleleDepths::type_info_gt_table[color].indexes[0]);
        haploid_prior(0) = haploid_prior_[color](pileup::AlleleDepths::type_info_table[color].indexes[0]);
        work_.SetFounders(diploid_prior, haploid_prior);
        work_.SetGenotypeLikelihoods(genotyper_, depths);
        double logdata = graph_.PeelForwards(work_, transition_matrices_[color]);
        prob_monomorphic_[color] = exp(logdata);
    }
}

// returns 'log10 P(Data ; model)-log10 scale' and log10 scaling.
LogProbability::value_t LogProbability::operator()(const pileup::RawDepths &depths,
                               int ref_index) {
    // calculate genotype likelihoods and store in the lower library vector
    double scale = work_.SetGenotypeLikelihoods(genotyper_, depths, ref_index);

    // Set the prior probability of the founders given the reference
    work_.SetFounders(diploid_prior_[ref_index], haploid_prior_[ref_index]);

    // for(int i=0; i < full_transition_matrices_.size(); ++i) {
    //     std::cerr << i << "\t";
    //     std::cerr << full_transition_matrices_[i].rows() << "x" << full_transition_matrices_[i].cols() << "\t";
    //     std::cerr << work_.upper[i].rows() << "x" << work_.upper[i].cols() << "\t";
    //     std::cerr << work_.lower[i].rows() << "x" << work_.lower[i].cols() << std::endl;
    // }

    // Calculate log P(Data ; model)
    double logdata = graph_.PeelForwards(work_, transition_matrices_[COLOR_ACGT]);

    return {logdata/M_LN10, scale/M_LN10};
}

// Calculate the probability of a depths object considering only indexes
// TODO: make indexes a property of the pedigree
LogProbability::value_t LogProbability::operator()(const pileup::AlleleDepths &depths) {
    const int ref_index = depths.type_info().reference;
    const int color = depths.color();

    double scale, logdata;
    // For monomorphic sites we have pre-calculated the peeling part
    if(color < 4) {
        assert(depths.type_info().width == 1 && depths.type_gt_info().width == 1 && ref_index == color);
        // Calculate genotype likelihoods
        scale = work_.SetGenotypeLikelihoods(genotyper_, depths);

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
        GenotypeArray diploid_prior = DiploidPrior(ref_index, color);
        GenotypeArray haploid_prior = HaploidPrior(ref_index, color);
        work_.SetFounders(diploid_prior, haploid_prior);
  
        // Calculate genotype likelihoods
        scale = work_.SetGenotypeLikelihoods(genotyper_, depths);

         // Calculate log P(Data ; model)
        logdata = graph_.PeelForwards(work_, transition_matrices_[color]);
    }
    return {logdata/M_LN10, scale/M_LN10};
}

// Construct the mutation matrices for each transition
TransitionMatrixVector dng::create_mutation_matrices(const RelationshipGraph &graph,
    const std::array<double, 4> &nuc_freq, const int mutype) {
    TransitionMatrixVector matrices(graph.num_nodes());
 
    for(size_t child = 0; child < graph.num_nodes(); ++child) {
        auto trans = graph.transition(child);
        if(trans.type == RelationshipGraph::TransitionType::Trio) {
            assert(graph.ploidy(child) == 2);
            auto dad = f81::matrix(trans.length1, nuc_freq);
            auto mom = f81::matrix(trans.length2, nuc_freq);
            matrices[child] = meiosis_matrix(graph.ploidy(trans.parent1), dad, graph.ploidy(trans.parent2), mom, mutype);
        } else if(trans.type == RelationshipGraph::TransitionType::Pair) {
            auto orig = f81::matrix(trans.length1, nuc_freq);
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

TransitionMatrix subset_mutation_matrix_meiosis_autosomal(const TransitionMatrix &full_matrix, size_t color) {
    const auto &type_info_gt_table = dng::pileup::AlleleDepths::type_info_gt_table;

    const int width = type_info_gt_table[color].width;

    TransitionMatrix ret{width*width, width};

    // a: parent1/dad; b: parent2/mom;
    // kronecker product is 'a-major'
    for(int a = 0; a < width; ++a) {
        const int ga = type_info_gt_table[color].indexes[a];
        for(int b = 0; b < width; ++b) {
            const int gb = type_info_gt_table[color].indexes[b];
            // copy correct value from the full matrix to the subset matrix
            const int x = a*width+b;
            const int gx = ga*10+gb;
            for(int y = 0; y < width; ++y) {
                const int gy = type_info_gt_table[color].indexes[y];
                // copy correct value from the full matrix to the subset matrix
                ret(x,y) = full_matrix(gx,gy); 
            }
        }
    }
    return ret;
}

TransitionMatrix subset_mutation_matrix_meiosis_xlinked(const TransitionMatrix &full_matrix, size_t color) {
    const auto &type_info_table = dng::pileup::AlleleDepths::type_info_table;
    const auto &type_info_gt_table = dng::pileup::AlleleDepths::type_info_gt_table;

    // a: parent1/dad; b: parent2/mom;
    // kronecker product is 'a-major'
    // a: haploid; b: diploid
    const int widthA = type_info_table[color].width;                      
    const int widthB = type_info_gt_table[color].width;

    TransitionMatrix ret{widthA*widthB,widthB};

    for(int a = 0; a < widthA; ++a) {
        const int ga = type_info_table[color].indexes[a];
        for(int b = 0; b < widthB; ++b) {
            const int gb = type_info_gt_table[color].indexes[b];
            // copy correct value from the full matrix to the subset matrix
            const int x = a*widthB+b;
            const int gx = ga*10+gb;
            for(int y = 0; y < widthB; ++y) {
                const int gy = type_info_gt_table[color].indexes[y];
                // copy correct value from the full matrix to the subset matrix
                ret(x,y) = full_matrix(gx,gy); 
            }
        }
    }
    return ret;
}

TransitionMatrix subset_mutation_matrix_meiosis_gamete(const TransitionMatrix &full_matrix, size_t color) {
    const auto &type_info_table = dng::pileup::AlleleDepths::type_info_table;
    const auto &type_info_gt_table = dng::pileup::AlleleDepths::type_info_gt_table;

    const int widthP = type_info_gt_table[color].width;
    const int widthC = type_info_table[color].width;

    TransitionMatrix ret{widthP,widthC};

    for(int x = 0; x < widthP; ++x) {
        const int gx = type_info_gt_table[color].indexes[x];
        for(int y = 0; y < widthC; ++y) {
            const int gy = type_info_table[color].indexes[y];
            // copy correct value from the full matrix to the subset matrix
            ret(x,y) = full_matrix(gx,gy); 
        }
    }
    return ret;
}


TransitionMatrix subset_mutation_matrix_mitosis_diploid(const TransitionMatrix &full_matrix, size_t color) {
    const auto &type_info_gt_table = dng::pileup::AlleleDepths::type_info_gt_table;

    const int width = type_info_gt_table[color].width;

    TransitionMatrix ret{width,width};

    for(int x = 0; x < width; ++x) {
        const int gx = type_info_gt_table[color].indexes[x];
        for(int y = 0; y < width; ++y) {
            const int gy = type_info_gt_table[color].indexes[y];
            // copy correct value from the full matrix to the subset matrix
            ret(x,y) = full_matrix(gx,gy); 
        }
    }           
    return ret;
}

TransitionMatrix subset_mutation_matrix_mitosis_haploid(const TransitionMatrix &full_matrix, size_t color) {
    const auto &type_info_table = dng::pileup::AlleleDepths::type_info_table;

    const int width = type_info_table[color].width;

    TransitionMatrix ret{width,width};

    for(int x = 0; x < width; ++x) {
        const int gx = type_info_table[color].indexes[x];
        for(int y = 0; y < width; ++y) {
            const int gy = type_info_table[color].indexes[y];
            // copy correct value from the full matrix to the subset matrix
            ret(x,y) = full_matrix(gx,gy); 
        }
    }
    return ret;
}


TransitionMatrixVector dng::create_mutation_matrices_subset(const TransitionMatrixVector &full_matrices, size_t color) {

    // Create out output vector
    TransitionMatrixVector matrices(full_matrices.size());

    // enumerate over all matrices
    for(size_t child = 0; child < full_matrices.size(); ++child) {
        if(full_matrices[child].rows() == 100 && full_matrices[child].cols() == 10) {
            // resize transition matrix to w*w,w
            matrices[child] = subset_mutation_matrix_meiosis_autosomal(full_matrices[child], color);
        } else if(full_matrices[child].rows() == 40 && full_matrices[child].cols() == 10) {
            // FIXME: ASSUME x-linkage for now
            matrices[child] = subset_mutation_matrix_meiosis_xlinked(full_matrices[child], color);
        } else if(full_matrices[child].rows() == 10 && full_matrices[child].cols() == 10) {
            matrices[child] = subset_mutation_matrix_mitosis_diploid(full_matrices[child], color);
        } else if(full_matrices[child].rows() == 4 && full_matrices[child].cols() == 4) {
            matrices[child] = subset_mutation_matrix_mitosis_haploid(full_matrices[child], color);
        } else if(full_matrices[child].rows() == 10 && full_matrices[child].cols() == 4) {
            matrices[child] = subset_mutation_matrix_meiosis_gamete(full_matrices[child], color);
        } else if(full_matrices[child].rows() == 0 && full_matrices[child].cols() == 0 ) {
            matrices[child] = {};
        } else {
            // should never reach here
            assert(false);
        }
    }
    return matrices;    
}
