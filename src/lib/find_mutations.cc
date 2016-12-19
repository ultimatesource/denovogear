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

FindMutations::FindMutations(double min_prob, const RelationshipGraph &graph,
        params_t params)
        : FindMutationsAbstract(min_prob, graph, params) {

    SetupPopulationPriorDiploid();
    SetupTransitionMatrix();
}

// Returns true if a mutation was found and the record was modified
bool FindMutations::operator()(const std::vector<depth_t> &depths,
        int ref_index, stats_t *stats) {
//    PR_NOTE(SW): Next version will be pass MutationStats instead of stats_t

    assert(stats != nullptr);
    //TODO(SW): Eventually, MutationStats this will replace all stats_t
    MutationStats mutation_stats(min_prob_);

    double scale = work_full_.SetGenotypeLikelihood(genotype_likelihood_,
                                                    depths, ref_index);

    work_full_.SetFounders(genotype_prior_[ref_index]);
    work_nomut_ = work_full_;

    bool is_mup_less_threshold = CalculateMutationProb(mutation_stats);

    if (is_mup_less_threshold) {
        return false;
    }

    relationship_graph_.PeelBackwards(work_full_, full_transition_matrices_);

    mutation_stats.SetGenotypeLikelihoods(work_full_, depths.size());
    mutation_stats.SetScaledLogLikelihood(scale);
    mutation_stats.SetPosteriorProbabilities(work_full_);

    mutation_stats.CalculateExpectedMutation(work_full_, mean_matrices_);
    mutation_stats.CalculateNodeMutation(work_full_,
                                         posmut_transition_matrices_);
    CalculateDenovoMutation(mutation_stats);

#if CALCULATE_ENTROPY == 1
    mutation_stats.CalculateEntropy(work_full_, posmut_transition_matrices_,
                                    max_entropies_, ref_index);
#endif

    //PR_NOTE(SW): Reassign mutation_stasts back to stats_t. This avoid major changes in call.cc at this stage
    //Remove this section after all PR are completed
    stats-> mup = mutation_stats.mup_;
    double pmut = stats->mup; //rest of the code will work
    stats-> lld = mutation_stats.lld_;
//    stats-> llh = mutation_stats.llh_;
    stats-> genotype_likelihoods = mutation_stats.genotype_likelihoods_;
    stats-> posterior_probabilities = mutation_stats.posterior_probabilities_;

	stats->mux = mutation_stats.mux_;

	stats->node_mup = mutation_stats.node_mup_;
	stats->mu1p = mutation_stats.mu1p_;

	stats->has_single_mut = mutation_stats.has_single_mut_;

	stats->dnq = mutation_stats.dnq_;
	stats->dnl = mutation_stats.dnl_;
	stats->dnt = mutation_stats.dnt_;

	stats->node_mu1p = mutation_stats.node_mu1p_;

#if CALCULATE_ENTROPY == 1
	stats->dnc = mutation_stats.dnc_;
#endif

    //End Remove this section

    return true;


}

void FindMutations::SetupTransitionMatrix(){

    for(size_t child = 0; child < work_full_.num_nodes; ++child) {
        auto trans = relationship_graph_.transitions()[child];

        if(trans.type == RelationshipGraph::TransitionType::Trio) {
            auto dad = f81::matrix(trans.length1, params_.nuc_freq);
            auto mom = f81::matrix(trans.length2, params_.nuc_freq);

            full_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom);
            nomut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 0);
            onemut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 1);
            mean_matrices_[child] = meiosis_diploid_mean_matrix(dad, mom);
        } else if(trans.type == RelationshipGraph::TransitionType::Pair) {
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);
            full_transition_matrices_[child] = mitosis_diploid_matrix(orig);
            nomut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 0);
            onemut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 1);
            mean_matrices_[child] = mitosis_diploid_mean_matrix(orig);
        }

        posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                             nomut_transition_matrices_[child];

    }

#if CALCULATE_ENTROPY == 1
    //Calculate max_entropy based on having no data
    for (int ref_index = 0; ref_index < 5; ++ref_index) {
        work_nomut_.SetFounders(genotype_prior_[ref_index]);

        relationship_graph_.PeelForwards(work_nomut_, nomut_transition_matrices_);
        relationship_graph_.PeelBackwards(work_nomut_, nomut_transition_matrices_);
        event_.assign(work_nomut_.num_nodes, 0.0);
        double total = 0.0, entropy = 0.0;
        for (std::size_t i = work_nomut_.founder_nodes.second;
                i < work_nomut_.num_nodes; ++i) {

            Eigen::ArrayXXd mat = (work_nomut_.super[i].matrix() *
                                   work_nomut_.lower[i].matrix().transpose()).array() *
                                   onemut_transition_matrices_[i].array();

            total += mat.sum();
            entropy += (mat.array() == 0.0).select(mat.array(),
                                                   mat.array() * mat.log()) .sum();
        }
        // Calculate entropy of mutation location
        max_entropies_[ref_index] = (-entropy / total + log(total)) / M_LN2;
    }
#endif

}
