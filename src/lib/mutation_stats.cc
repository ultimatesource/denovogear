/*
 * Copyright (c) 2014-2016 Reed A. Cartwright
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *           Steven H. Wu <stevenwu@asu.edu>
 *
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

#include <dng/mutation_stats.h>

MutationStats::MutationStats(double min_prob)
        : min_prob_(min_prob) {
}

bool MutationStats::CalculateMutationProb(
        const dng::peel::workspace_t &work_nomut,
        const dng::peel::workspace_t &work_full) {
	// P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
	logdata_nomut_ = work_nomut.forward_result;
	logdata_ = work_full.forward_result;
	mup_ = static_cast<float>(-std::expm1(logdata_nomut_ - logdata_));
	return mup_ < min_prob_;
}

void MutationStats::SetScaledLogLikelihood(double scale) {
	lld_ = static_cast<float>((logdata_ + scale) / M_LN10);
//    llh_ = static_cast<float>( logdata_ / M_LN10);
}

void MutationStats::SetGenotypeLikelihoods(
		const dng::peel::workspace_t &workspace, const int depth_size) {

	genotype_likelihoods_.resize(workspace.num_nodes);
	for (std::size_t u = 0; u < depth_size; ++u) {
		std::size_t pos = workspace.library_nodes.first + u;
		genotype_likelihoods_[pos] = workspace.lower[pos].log() / M_LN10;
		//TODO(SW): Eigen 3.3 might have log10()
	}
}

void MutationStats::SetPosteriorProbabilities(
		const dng::peel::workspace_t &workspace) {

	posterior_probabilities_.resize(workspace.num_nodes);
	for (std::size_t i = 0; i < workspace.num_nodes; ++i) {
//        posterior_probabilities_[i] = WORKSPACE_MULTIPLE_UPPER_LOWER(workspace, i);
		posterior_probabilities_[i] = (workspace.upper[i] * workspace.lower[i]);
		posterior_probabilities_[i] /= posterior_probabilities_[i].sum();
	}

}




void MutationStats::CalculateExpectedMutation(dng::peel::workspace_t &work_full,
                                              dng::TransitionVector &mean_matrices){
    mux_ = 0.0;
    for(size_t i = work_full.founder_nodes.second; i < work_full.num_nodes; ++i) {
        mux_ += (work_full.super[i] * (mean_matrices[i] *
                                      work_full.lower[i].matrix()).array()).sum();
    }

};

void MutationStats::CalculateNodeMutation(dng::peel::workspace_t &work_full,
                                          dng::TransitionVector &posmut_transition_matrices) {
    std::vector<double> event (work_full.num_nodes, 0.0);
    for (size_t i = work_full.founder_nodes.second; i < work_full.num_nodes; ++i) {
        event[i] = (work_full.super[i] * (posmut_transition_matrices[i] *
                                          work_full.lower[i].matrix()).array()).sum();
        event[i] = event[i] / mup_;
    }

    SetNodeMup(event, work_full.founder_nodes.second);


}

void MutationStats::CalculateDenovoMutation(dng::peel::workspace_t &work_nomut,
                                            dng::TransitionVector &onemut_transition_matrices,
                                            const dng::RelationshipGraph &pedigree) {
    std::vector<double> event (work_nomut.num_nodes, 0.0);
    double total = 0.0, max_coeff = -1.0;
    size_t dn_row = 0, dn_col = 0, dn_location = 0;

    for (std::size_t i = work_nomut.founder_nodes.second; i < work_nomut.num_nodes; ++i) {
        Eigen::ArrayXXd mat = (work_nomut.super[i].matrix() *
                               work_nomut.lower[i].matrix().transpose()).array() *
                              onemut_transition_matrices[i].array();

        UpdateMaxDeNovoMutation(mat, i, max_coeff, dn_row, dn_col, dn_location);

        event[i] = mat.sum();
        total += event[i];

    }
    SetExactlyOneMutation(total);

    if (has_single_mut_) {
        SetNodeMu1p(event, total, work_nomut.founder_nodes.second);

        dnq_ = dng::utility::lphred<int32_t>(1.0 - (max_coeff / total), 255);
        dnl_ = pedigree.labels()[dn_location];
        if (pedigree.transitions()[dn_location].type == dng::RelationshipGraph::TransitionType::Germline) {
            dnt_ = &dng::meiotic_diploid_mutation_labels[dn_row][dn_col][0];
        } else {
            dnt_ = &dng::mitotic_diploid_mutation_labels[dn_row][dn_col][0];
        }

    }
}


#if CALCULATE_ENTROPY == 1
void MutationStats::CalculateEntropy(dng::peel::workspace_t &work_nomut,
                                     dng::TransitionVector &onemut_transition_matrices,
                                     std::array<double, 5> max_entropies,
                                     std::size_t ref_index) {

    double total = 0.0;
    double entropy = 0.0;

    if (has_single_mut_) {

        for (std::size_t i = work_nomut.founder_nodes.second; i < work_nomut.num_nodes; ++i) {
            Eigen::ArrayXXd mat = (work_nomut.super[i].matrix() *
                                   work_nomut.lower[i].matrix().transpose()).array() *
                                  onemut_transition_matrices[i].array();
            total += mat.sum();

            entropy += (mat.array() == 0.0).select(mat.array(),
                                                   mat.array() * mat.log()).sum();
        }

        entropy = (-entropy / total + log(total)) / M_LN2;
        entropy /= max_entropies[ref_index];
        dnc_ = static_cast<float>(std::round(100.0 * (1.0 - entropy)));

    }

}
#endif

void MutationStats::SetNodeMup(const std::vector<double> &event,
                               std::size_t first_nonfounder_index) {
    SetNodeCore(node_mup_, event, first_nonfounder_index);
}

void MutationStats::SetNodeMu1p(std::vector<double> &event, double total,
                                std::size_t first_nonfounder_index) {
    for (std::size_t i = first_nonfounder_index; i < event.size(); ++i) {
        event[i] = event[i] / total;
    }
    SetNodeCore(node_mu1p_, event, first_nonfounder_index);
}

void MutationStats::SetNodeCore(std::vector<float> &stats,
                                const std::vector<double> &event,
                                std::size_t first_nonfounder_index) {
    stats.resize(event.size(), hts::bcf::float_missing);
    for (std::size_t i = first_nonfounder_index; i < event.size(); ++i) {
        stats[i] = static_cast<float>(event[i]);
    }
    //HOWTO: use std::copy with cast??    std::copy(event.begin()+first_nonfounder_index, event.end(), stats.begin() );
}



void MutationStats::SetExactlyOneMutation(double total){
    mu1p_ = static_cast<float>(total * (1.0 - mup_));
    has_single_mut_ = ((mu1p_ / mup_) >= min_prob_);
}

void MutationStats::UpdateMaxDeNovoMutation(const Eigen::ArrayXXd &mat,
                                            size_t &index, double &max_coeff,
                                            size_t &dn_row, size_t &dn_col,
                                            size_t &dn_location) {

    std::size_t row, col;
    double mat_max = mat.maxCoeff(&row, &col);
    if (mat_max > max_coeff) {
        max_coeff = mat_max;
        dn_row = row;
        dn_col = col;
        dn_location = index;
    }

}
