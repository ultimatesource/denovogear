//
// Created by steven on 1/15/16.
//

#ifndef DENOVOGEAR_FIND_MUTATION_GETTER_H
#define DENOVOGEAR_FIND_MUTATION_GETTER_H


#include "find_mutation.h"

//Just make testing easier, no real purpose of this class
class FindMutationsGetter : public FindMutations {


public:
    FindMutationsGetter(double min_prob, const Pedigree &pedigree, const params_t &params) : FindMutations(min_prob,
                                                                                                           pedigree,
                                                                                                           params) { }
public:

    const Pedigree &getPedigree_() const {
        return pedigree_;
    }

    const params_t &getParams_() const {
        return params_;
    }

    double getMin_prob_() const {
        return min_prob_;
    }

    const peel::workspace_t &getworkspace() const {
        return work_;
    }

    const TransitionVector &getFull_transition_matrices_() const {
        return full_transition_matrices_;
    }

    const TransitionVector &getNomut_transition_matrices_() const {
        return nomut_transition_matrices_;
    }

    const TransitionVector &getPosmut_transition_matrices_() const {
        return posmut_transition_matrices_;
    }

    const TransitionVector &getOnemut_transition_matrices_() const {
        return onemut_transition_matrices_;
    }

    const TransitionVector &getMean_matrices_() const {
        return mean_matrices_;
    }

    const genotype::DirichletMultinomialMixture &getGenotype_likelihood_() const {
        return genotype_likelihood_;
    }

    const GenotypeArray *getGenotype_prior_() const {
        return genotype_prior_;
    }

    const std::array<double, 5> &getMax_entropies_() const {
        return max_entropies_;
    }

    const std::vector<double> &getEvent_() const {
        return event_;
    }

};

#endif //DENOVOGEAR_FIND_MUTATION_GETTER_H
