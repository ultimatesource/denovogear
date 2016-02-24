//
// Created by steven on 2/3/16.
//

#ifndef DENOVOGEAR_MUTATION_STATS_H
#define DENOVOGEAR_MUTATION_STATS_H

#include <string>
#include <vector>

#include <dng/matrix.h>
#include <dng/peeling.h>
#include <dng/hts/bcf.h>

class MutationStats {

public:

    MutationStats(double min_prob);

    float get_mutation_prob() const;

    bool get_has_single_mut() const;

    const dng::GenotypeArray &inspect_posterior_at(int index) const;

    const dng::GenotypeArray &inspect_genotype_at(int index) const;


    bool set_mutation_prob(const double logdata_nomut, const double logdata);

    void set_scaled_log_likelihood(double scale);

    void set_genotype_likelihood(const dng::peel::workspace_t &workspace, const int depth_size);

    void set_posterior_probabilities(const dng::peel::workspace_t &workspace);

    void set_exactly_one_mutation(double total);

    void set_node_mup(const std::vector<double> &event, std::size_t first_nonfounder_index);

    void set_node_mu1p(std::vector<double> &event, double total, std::size_t first_nonfounder_index);



    void record_basic_stats(hts::bcf::Variant &record);

    void record_genotype_stats(hts::bcf::Variant &record);

    void record_single_mutation_stats(hts::bcf::Variant &record);


    enum class OutputLevel {BASIC, COMPLETE, SINGLE};
    //TODO: Record different sets of variables, [basic, complete, everytihng/single], something like the following
    //basic [ mup, lld, llh, genotype_likelihood]?
    //complete [basic, mux, posterior, node_mup]
    //single [complete, has_single_mut, mu1p, dnt, dnl, dnq, dnc, node_mu1p]



    void set_genotype_related_stats(const int (&acgt_to_refalt_allele)[5],
                                                   const int (&refalt_to_acgt_allele)[5],
                                                   const uint32_t n_alleles,
                                                   const std::size_t ref_index,
                                                   const std::size_t num_nodes,
                                                   const std::size_t library_start);



private://TODO: surely these should not be public. Refactor these while working with call.cc
public:
    float mup;
    float lld;
    float llh;
    float mux;

    bool has_single_mut;
    float mu1p;

    std::string dnt;
    std::string dnl;
    int32_t dnq;
    [[deprecated]] int32_t dnc;

    dng::IndividualVector posterior_probabilities;
    dng::IndividualVector genotype_likelihoods;
    std::vector<float> node_mup;
    std::vector<float> node_mu1p;

    std::vector<int32_t> best_genotypes;//.resize(2 * num_nodes);
    std::vector<int32_t> genotype_qualities;//(num_nodes);
    std::vector<float> gp_scores;//(gt_count*num_nodes );
    std::vector<float> gl_scores;//(gt_count *num_nodes, hts::bcf::float_missing);

private:

    double min_prob;
    double logdata;//TODO: Maybe hack these away, using lld, llh
    double logdata_nomut;

    void set_node_core(std::vector<float> &stats, const std::vector<double> &event,
                       std::size_t first_nonfounder_index);

};
#endif //DENOVOGEAR_MUTATION_STATS_H
