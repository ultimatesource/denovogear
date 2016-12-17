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


#include <dng/find_mutations_ylinked.h>


using namespace dng;

FindMutationsYLinked::~FindMutationsYLinked() {
	// TODO Auto-generated destructor stub
}


FindMutationsYLinked::FindMutationsYLinked(double min_prob,
        const RelationshipGraph &graph, params_t params)
        : FindMutationsAbstract(min_prob, graph, params)
    {

    SetupPopulationPriorHaploid();
    SetupTransitionMatrix();

#if CALCULATE_ENTROPY == 1
    std::cerr << "Entropy will not be calculated for Y-linked model!!" << std::endl;
#endif

}

// Returns true if a mutation was found and the record was modified
bool FindMutationsYLinked::operator()(const std::vector<depth_t> &depths,
                               int ref_index, stats_t *stats) {
    using namespace hts::bcf;
    using dng::utility::lphred;
    using dng::utility::phred;

    assert(stats != nullptr);

    //TODO(SW): Eventually, MutationStats this will replace all stats_t
    MutationStats mutation_stats(min_prob_);

    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::tie(work_full_.lower[work_full_.library_nodes.first + u], stemp) =
                genotype_likelihood_.CalculateHaploid(depths[u], ref_index);
        scale += stemp;
    }
    // Set the prior probability of the founders given the reference
    work_full_.SetFounders(genotype_prior_[ref_index]);
    work_nomut_ = work_full_;

    bool is_mup_less_threshold = CalculateMutationProb(mutation_stats);

    //TODO(SW): HACK: for unittest
    stats-> mup = mutation_stats.mup_;
    if (is_mup_less_threshold) {
        return false;
    }

    relationship_graph_.PeelBackwards(work_full_, full_transition_matrices_);

    mutation_stats.SetGenotypeLikelihoods(work_full_, depths.size());
    mutation_stats.SetScaledLogLikelihood(scale);
    mutation_stats.SetPosteriorProbabilities(work_full_);

    //PR_NOTE(SW): Reassign mutation_stasts back to stats_t. This avoid major changes in call.cc at this stage
    //Remove this section after all PR are completed
    stats-> mup = mutation_stats.mup_;
    double pmut = stats->mup; //rest of the code will work
    stats-> lld = mutation_stats.lld_;
//    stats-> llh = mutation_stats.llh_;
    stats-> genotype_likelihoods = mutation_stats.genotype_likelihoods_;
    stats-> posterior_probabilities = mutation_stats.posterior_probabilities_;
    //End Remove this section

    //TODO(SW): Haven't test the following
    mutation_stats.CalculateExpectedMutation(work_full_, mean_matrices_);
    mutation_stats.CalculateNodeMutation(work_full_,
                                         posmut_transition_matrices_);
    CalculateDenovoMutation(mutation_stats);

    stats->mux = mutation_stats.mux_;
    stats->node_mup = mutation_stats.node_mup_;
    stats->mu1p = mutation_stats.mu1p_;
    stats->has_single_mut = mutation_stats.has_single_mut_;
    stats->dnq = mutation_stats.dnq_;
    stats->dnl = mutation_stats.dnl_;
    stats->dnt = mutation_stats.dnt_;
    stats->node_mu1p = mutation_stats.node_mu1p_;


    return true;
}



void FindMutationsYLinked::SetupTransitionMatrix(){

    for(size_t child = 0; child < work_full_.num_nodes; ++child) {
        auto trans = relationship_graph_.transitions()[child];

        if (trans.type == RelationshipGraph::TransitionType::Pair) {
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);
            full_transition_matrices_[child] = mitosis_haploid_matrix(orig);
            nomut_transition_matrices_[child] = mitosis_haploid_matrix(orig, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = mitosis_haploid_matrix(orig, 1);
            mean_matrices_[child] = mitosis_haploid_mean_matrix(orig);
        }
    }


}
