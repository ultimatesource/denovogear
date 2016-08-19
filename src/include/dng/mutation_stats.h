/*
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
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
#ifndef DENOVOGEAR_MUTATION_STATS_H
#define DENOVOGEAR_MUTATION_STATS_H

#include <string>
#include <vector>

#include <dng/hts/bcf.h>
#include <dng/matrix.h>
#include <dng/peeling.h>
#include <dng/pedigree.h>
#include <dng/mutation.h>
#include <dng/detail/unit_test.h>

//TODO(SW): Need to change all test cases before remove entropy calculation
#define CALCULATE_ENTROPY 1

class MutationStats {

public:

	MutationStats(double min_prob);

	bool CalculateMutationProb(const dng::peel::workspace_t &work_nomut,
			const dng::peel::workspace_t &work_full);

	void SetScaledLogLikelihood(double scale);

	void SetGenotypeLikelihoods(const dng::peel::workspace_t &workspace,
			const int depth_size);

	void SetPosteriorProbabilities(const dng::peel::workspace_t &workspace);

//PR_NOTE(SW): Functions for the next PR
//	void CalculateExpectedMutation(dng::peel::workspace_t &work_full,
//			dng::TransitionVector &mean_matrices);
//
//	void CalculateNodeMutation(dng::peel::workspace_t &work_full,
//			dng::TransitionVector &posmut_transition_matrices);
//
//	void CalculateDenovoMutation(dng::peel::workspace_t &work_nomut,
//			dng::TransitionVector &onemut_transition_matrices,
//			const dng::Pedigree &pedigree);

public:
	//TODO(SW): think about whether these should be public or private?
	float mup_;
	float lld_;
	[[deprecated]] float llh_;
	float lls_;
	float mux_;

	bool has_single_mut_;
	float mu1p_;

	std::string dnt_;
	std::string dnl_;
	int32_t dnq_;
	[[deprecated]]int32_t dnc_;

	dng::IndividualVector posterior_probabilities_;
	dng::IndividualVector genotype_likelihoods_;
	std::vector<float> node_mup_;
	std::vector<float> node_mu1p_;

	std::vector<int32_t> best_genotypes_; //.resize(2 * num_nodes);
	std::vector<int32_t> genotype_qualities_; //(num_nodes);
	std::vector<float> gp_scores_; //(gt_count * num_nodes);
	std::vector<float> gl_scores_; //(gt_count * num_nodes, hts::bcf::float_missing);

private:

	double min_prob_;
	double logdata_;
	double logdata_nomut_;

//PR_NOTE(SW): Functions for the next PR
//	void SetExactlyOneMutation(double total);
//
//	void SetNodeCore(std::vector<float> &stats,
//			const std::vector<double> &event,
//			std::size_t first_nonfounder_index);
//
//	void UpdateMaxDeNovoMutation(const Eigen::ArrayXXd &mat, size_t &index,
//			double &max_coeff, size_t &dn_row, size_t &dn_col, size_t &dn_loc);
//#if CALCULATE_ENTROPY == 1
//public:
//	[[deprecated]]
//	void CalculateEntropy(dng::peel::workspace_t &work_nomut,
//			dng::TransitionVector &onemut_transition_matrices,
//			std::array<double, 5> max_entropies, std::size_t ref_index);
//#endif

};

#endif //DENOVOGEAR_MUTATION_STATS_H
