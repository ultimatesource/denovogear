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


#include <dng/find_mutations.h>

using namespace dng;

FindMutations::~FindMutations() {
	// TODO Auto-generated destructor stub
}


FindMutations::FindMutations(double min_prob, const Pedigree &pedigree,
                             params_t params) :
    pedigree_(pedigree), min_prob_(min_prob),
    params_(params), genotype_likelihood_{params.params_a, params.params_b},
    work_(pedigree.CreateWorkspace()) {

    using namespace dng;

    // Use a parent-independent mutation model, which produces a
    // beta-binomial
    genotype_prior_[0] = population_prior(params_.theta, params_.nuc_freq, {params_.ref_weight, 0, 0, 0});
    genotype_prior_[1] = population_prior(params_.theta, params_.nuc_freq, {0, params_.ref_weight, 0, 0});
    genotype_prior_[2] = population_prior(params_.theta, params_.nuc_freq, {0, 0, params_.ref_weight, 0});
    genotype_prior_[3] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, params_.ref_weight});
    genotype_prior_[4] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, 0});

    // Calculate mutation expectation matrices
    full_transition_matrices_.assign(work_.num_nodes, {});
    nomut_transition_matrices_.assign(work_.num_nodes, {});
    posmut_transition_matrices_.assign(work_.num_nodes, {});
    onemut_transition_matrices_.assign(work_.num_nodes, {});
    mean_matrices_.assign(work_.num_nodes, {});

    for(size_t child = 0; child < work_.num_nodes; ++child) {
        auto trans = pedigree.transitions()[child];
        if(trans.type == Pedigree::TransitionType::Germline) {
            auto dad = f81::matrix(trans.length1, params_.nuc_freq);
            auto mom = f81::matrix(trans.length2, params_.nuc_freq);

            full_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom);
            nomut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 1);
            mean_matrices_[child] = meiosis_diploid_mean_matrix(dad, mom);
        } else if(trans.type == Pedigree::TransitionType::Somatic ||
                  trans.type == Pedigree::TransitionType::Library) {
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);

            full_transition_matrices_[child] = mitosis_diploid_matrix(orig);
            nomut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 1);
            mean_matrices_[child] = mitosis_diploid_mean_matrix(orig);
        } else {
            full_transition_matrices_[child] = {};
            nomut_transition_matrices_[child] = {};
            posmut_transition_matrices_[child] = {};
            onemut_transition_matrices_[child] = {};
            mean_matrices_[child] = {};
        }
    }

    //Calculate max_entropy based on having no data
    for(int ref_index = 0; ref_index < 5; ++ref_index) {
        work_.SetFounders(genotype_prior_[ref_index]);

        pedigree_.PeelForwards(work_, nomut_transition_matrices_);
        pedigree_.PeelBackwards(work_, nomut_transition_matrices_);
        event_.assign(work_.num_nodes, 0.0);
        double total = 0.0, entropy = 0.0;
        for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
            Eigen::ArrayXXd mat = (work_.super[i].matrix() *
                                   work_.lower[i].matrix().transpose()).array() *
                                  onemut_transition_matrices_[i].array();
            total += mat.sum();
            entropy += (mat.array() == 0.0).select(mat.array(),
                                                   mat.array() * mat.log()).sum();
        }
        // Calculate entropy of mutation location
        max_entropies_[ref_index] = (-entropy / total + log(total)) / M_LN2;
    }
}

// Returns true if a mutation was found and the record was modified
bool FindMutations::operator()(const std::vector<depth_t> &depths,
                               int ref_index, stats_t *stats) {
    using namespace std;
    using namespace hts::bcf;
    using dng::utility::lphred;
    using dng::utility::phred;

    assert(stats != nullptr);

    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::tie(work_.lower[work_.library_nodes.first + u], stemp) =
            genotype_likelihood_(depths[u], ref_index);
        scale += stemp;
    }

    // Set the prior probability of the founders given the reference
    work_.SetFounders(genotype_prior_[ref_index]);

    // Calculate log P(Data, nomut ; model)
    const double logdata_nomut = pedigree_.PeelForwards(work_,
                                 nomut_transition_matrices_);

    /**** Forward-Backwards with full-mutation ****/

    // Calculate log P(Data ; model)
    const double logdata = pedigree_.PeelForwards(work_, full_transition_matrices_);

    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    const double pmut = -std::expm1(logdata_nomut - logdata);

    // Skip this site if it does not meet lower probability threshold
    if(pmut < min_prob_) {
        return false;
    }

    // Copy genotype likelihoods
    stats->genotype_likelihoods.resize(work_.num_nodes);
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::size_t pos = work_.library_nodes.first + u;
        stats->genotype_likelihoods[pos] = work_.lower[pos].log() / M_LN10;
    }

    // Peel Backwards with full-mutation
    pedigree_.PeelBackwards(work_, full_transition_matrices_);

    size_t library_start = work_.library_nodes.first;

    stats->mup = pmut;
    stats->lld = (logdata + scale) / M_LN10;
//    stats->llh = logdata / M_LN10;

    // Calculate statistics after Forward-Backwards
    stats->posterior_probabilities.resize(work_.num_nodes);
    for(std::size_t i = 0; i < work_.num_nodes; ++i) {
        stats->posterior_probabilities[i] = work_.upper[i] * work_.lower[i];
        stats->posterior_probabilities[i] /= stats->posterior_probabilities[i].sum();
    }
    double mux = 0.0;
    event_.assign(work_.num_nodes, 0.0);
    for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        mux += (work_.super[i] * (mean_matrices_[i] *
                                  work_.lower[i].matrix()).array()).sum();
        event_[i] = (work_.super[i] * (posmut_transition_matrices_[i] *
                                       work_.lower[i].matrix()).array()).sum();
        event_[i] = event_[i] / pmut;
    }
    stats->mux = mux;

    stats->node_mup.resize(work_.num_nodes, hts::bcf::float_missing);
    for(size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        stats->node_mup[i] = static_cast<float>(event_[i]);
    }

    /**** Forward-Backwards with no-mutation ****/

    // TODO: Better to use a separate workspace???
    pedigree_.PeelForwards(work_, nomut_transition_matrices_);
    pedigree_.PeelBackwards(work_, nomut_transition_matrices_);
    event_.assign(work_.num_nodes, 0.0);
    double total = 0.0, entropy = 0.0, max_coeff = -1.0;
    size_t dn_loc = 0, dn_col = 0, dn_row = 0;
    for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        Eigen::ArrayXXd mat = (work_.super[i].matrix() *
                               work_.lower[i].matrix().transpose()).array() *
                              onemut_transition_matrices_[i].array();
        std::size_t row, col;
        double mat_max = mat.maxCoeff(&row, &col);
        if(mat_max > max_coeff) {
            max_coeff = mat_max;
            dn_row  = row;
            dn_col = col;
            dn_loc = i;
        }
        event_[i] = mat.sum();
        entropy += (mat.array() == 0.0).select(mat.array(),
                                               mat.array() * mat.log()).sum();
        total += event_[i];
    }
    // Calculate P(only 1 mutation)
    const double pmut1 = total * (1.0 - pmut);
    stats->mu1p = pmut1;

    // Output statistics for single mutation only if it is likely
    if(pmut1 / pmut >= min_prob_) {
        stats->has_single_mut = true;
        for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
            event_[i] = event_[i] / total;
        }

        // Calculate entropy of mutation location
        entropy = (-entropy / total + log(total)) / M_LN2;
        entropy /= max_entropies_[ref_index];
        stats->dnc = std::round(100.0 * (1.0 - entropy));

        stats->dnq = lphred<int32_t>(1.0 - (max_coeff / total), 255);
        stats->dnl = pedigree_.labels()[dn_loc];
        if(pedigree_.transitions()[dn_loc].type == Pedigree::TransitionType::Germline) {
            stats->dnt = &meiotic_diploid_mutation_labels[dn_row][dn_col][0];
        } else {
            stats->dnt = &mitotic_diploid_mutation_labels[dn_row][dn_col][0];
        }

        stats->node_mu1p.resize(work_.num_nodes, hts::bcf::float_missing);
        for(size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
            stats->node_mu1p[i] = static_cast<float>(event_[i]);
        }
    } else {
        stats->has_single_mut = false;
    }
    return true;
}
