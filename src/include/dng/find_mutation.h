//
// Created by steven on 1/11/16.
//

#pragma once
#ifndef DENOVOGEAR_FIND_MUTATION_H
#define DENOVOGEAR_FIND_MUTATION_H

#include <iostream>

#include <dng/hts/bam.h>
#include <dng/hts/bcf.h>
#include <dng/pedigree.h>
#include <dng/likelihood.h>
#include <dng/mutation.h>

#include "mutation_stats.h"

//using namespace dng::task;
using namespace dng;

std::vector<std::pair<std::string, uint32_t>> parse_contigs(const bam_hdr_t *hdr);

std::vector<std::string> extract_contigs(const bcf_hdr_t *hdr) ;

class FindMutations {
public:
    struct params_t {
        double theta;
        std::array<double, 4> nuc_freq;
        double ref_weight;

        dng::genotype::DirichletMultinomialMixture::params_t params_a;
        dng::genotype::DirichletMultinomialMixture::params_t params_b;
    };

    //TODO: replace struct stats_t with mutation_stats.cc
//    struct [[deprecated]] stats_t {
    struct stats_t {
        float mup;
        float lld;
        float llh;
        float mux;

        bool has_single_mut;
        float mu1p;
        std::string dnt;
        std::string dnl;
        int32_t dnq;
        int32_t dnc;

        IndividualVector posterior_probabilities;
        IndividualVector genotype_likelihoods;
        std::vector<float> node_mup;
        std::vector<float> node_mu1p;
    };

    FindMutations(double min_prob, const Pedigree &pedigree, params_t params);

    bool operator()(const std::vector<depth_t> &depths, int ref_index,
                    stats_t *stats);

    //TODO: either place with this function, or replace operator() with this
    bool calculate_mutation(const std::vector<depth_t> &depths, int ref_index, MutationStats &mutation_stats);

protected:

    const dng::Pedigree &pedigree_;
    params_t params_;
    double min_prob_;
    dng::peel::workspace_t work_;

    dng::TransitionVector full_transition_matrices_;
    dng::TransitionVector nomut_transition_matrices_;
    dng::TransitionVector posmut_transition_matrices_;
    dng::TransitionVector onemut_transition_matrices_;
    dng::TransitionVector mean_matrices_;

    // Model genotype likelihoods as a mixture of two dirichlet multinomials
    // TODO: control these with parameters
    dng::genotype::DirichletMultinomialMixture genotype_likelihood_;

    dng::GenotypeArray genotype_prior_[5]; // Holds P(G | theta)
    std::array<double, 5> max_entropies_;
    std::vector<double> event_;


    bool calculate_mutation_prob(MutationStats &stats);
    void calculate_posterior_probabilities(MutationStats &mutation_stats);
    void calculate_expected_mutation(MutationStats &mutation_stats);
    void calculate_node_mutation(MutationStats &mutation_stats);
    void calculate_denovo_mutation(MutationStats &mutation_stats);


};

#endif //DENOVOGEAR_FIND_MUTATION_H
