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


#include <dng/find_mutations_abstract.h>


using namespace dng;
//TODO(SW): Eventually this will be just FindMutationsAbstract.
//TODO(SW): The original FindMutationsAbstract will become FindMutationsAbstractAutosomal

FindMutationsAbstract::~FindMutationsAbstract() {
	// TODO Auto-generated destructor stub
}


FindMutationsAbstract::FindMutationsAbstract(double min_prob, const RelationshipGraph &graph,
                             params_t params) :
    relationship_graph_(graph), min_prob_(min_prob),
    params_(params), genotype_likelihood_{params.params_a, params.params_b},
    work_full_(graph.CreateWorkspace()), work_nomut_(graph.CreateWorkspace()) {

    // Calculate mutation expectation matrices
    full_transition_matrices_.assign(work_full_.num_nodes, {});
    nomut_transition_matrices_.assign(work_full_.num_nodes, {});
    posmut_transition_matrices_.assign(work_full_.num_nodes, {});
    onemut_transition_matrices_.assign(work_full_.num_nodes, {});
    mean_matrices_.assign(work_full_.num_nodes, {});

    //keep_library_index_ = graph.KeepLibraryIndex();
}

void FindMutationsAbstract::SetupPopulationPriorDiploid() {
    genotype_prior_[0] = population_prior_diploid(params_.theta, params_.nuc_freq,
                                          {params_.ref_weight, 0, 0, 0});
    genotype_prior_[1] = population_prior_diploid(params_.theta, params_.nuc_freq,
                                          {0, params_.ref_weight, 0, 0});
    genotype_prior_[2] = population_prior_diploid(params_.theta, params_.nuc_freq,
                                          {0, 0, params_.ref_weight, 0});
    genotype_prior_[3] = population_prior_diploid(params_.theta, params_.nuc_freq,
                                          {0, 0, 0, params_.ref_weight});
    genotype_prior_[4] = population_prior_diploid(params_.theta, params_.nuc_freq,
                                          {0, 0, 0, 0});
}

void FindMutationsAbstract::SetupPopulationPriorHaploid() {
    genotype_prior_[0] = population_prior_haploid(params_.theta, params_.nuc_freq,
                                                {params_.ref_weight, 0, 0, 0});
    genotype_prior_[1] = population_prior_haploid(params_.theta, params_.nuc_freq,
                                                {0, params_.ref_weight, 0, 0});
    genotype_prior_[2] = population_prior_haploid(params_.theta, params_.nuc_freq,
                                                {0, 0, params_.ref_weight, 0});
    genotype_prior_[3] = population_prior_haploid(params_.theta, params_.nuc_freq,
                                                {0, 0, 0, params_.ref_weight});
    genotype_prior_[4] = population_prior_haploid(params_.theta, params_.nuc_freq,
                                                {0, 0, 0, 0});
}

void FindMutationsAbstract::Resize10To4(GenotypeArray &array){

    array(1, 0) = array(4, 0);
    array(2, 0) = array(7, 0);
    array(3, 0) = array(9, 0);

    array = array.block<4,1>(0,0).eval();


}

bool FindMutationsAbstract::CalculateMutationProb(MutationStats &mutation_stats) {

    // Calculate log P(Data, nomut ; model)
	relationship_graph_.PeelForwards(work_nomut_, nomut_transition_matrices_);

	// Calculate log P(Data ; model)
	relationship_graph_.PeelForwards(work_full_, full_transition_matrices_);
	// P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]

	bool is_mup_less_threshold = mutation_stats.CalculateMutationProb(
			work_nomut_, work_full_);

	return is_mup_less_threshold;
}


void FindMutationsAbstract::CalculateDenovoMutation(MutationStats &mutation_stats) {

    relationship_graph_.PeelBackwards(work_nomut_, nomut_transition_matrices_);
    mutation_stats.CalculateDenovoMutation(work_nomut_,
                                           onemut_transition_matrices_,
                                           relationship_graph_);

}






