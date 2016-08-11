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

MutationStats::MutationStats(double min_prob) :
		min_prob_(min_prob) {
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
