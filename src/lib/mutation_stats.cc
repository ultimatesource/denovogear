//
// Created by steven on 2/3/16.
//
#include <iostream>

#include <dng/peeling.h>

#include <dng/mutation_stats.h>


MutationStats::MutationStats(double min_prob) : min_prob(min_prob) { }


bool MutationStats::set_mutation_prob(const double logdata_nomut, const double logdata){
    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    this->logdata_nomut = logdata_nomut;
    this->logdata = logdata;
    mup = -std::expm1(logdata_nomut - logdata) ;
    return mup < min_prob;
}


void MutationStats::set_scaled_log_likelihood(double scale) {
    lld = (logdata + scale) / M_LN10;
    llh = logdata / M_LN10;

}

void MutationStats::set_genotype_likelihood(const dng::peel::workspace_t &workspace, const int depth_size) {

    genotype_likelihoods.resize(workspace.num_nodes);
    for(std::size_t u = 0; u < depth_size; ++u) {
        std::size_t pos = workspace.library_nodes.first + u;
        genotype_likelihoods[pos] = workspace.lower[pos].log() / M_LN10;
        //TODO: Eigen 3.3 might have log10()
    }
}

void MutationStats::set_posterior_probabilities(const dng::peel::workspace_t &workspace) {

    posterior_probabilities.resize(workspace.num_nodes);
    for(std::size_t i = 0; i < workspace.num_nodes; ++i) {
        posterior_probabilities[i] = WORKSPACE_T_MULTIPLE_UPPER_LOWER(workspace, i);
        posterior_probabilities[i] /= posterior_probabilities[i].sum();
    }

}

void MutationStats::set_exactly_one_mutation(double total){
    mu1p = total * (1.0 - mup);
    has_single_mut = (mu1p / mup) >= min_prob;
}

void MutationStats::set_node_mup(const std::vector<double> &event,
                                 std::size_t first_nonfounder_index) {
    set_node_core(node_mup, event, first_nonfounder_index);
}

void MutationStats::set_node_mu1p(std::vector<double> &event, double total,
                                  std::size_t first_nonfounder_index) {
    for (std::size_t i = first_nonfounder_index; i < event.size(); ++i) {
        event[i] = event[i] / total;
    }
    set_node_core(node_mu1p, event, first_nonfounder_index);
}

void MutationStats::set_node_core(std::vector<float> &stats, const std::vector<double> &event,
                                  std::size_t first_nonfounder_index) {
    stats.resize(event.size(), hts::bcf::float_missing);
    for (std::size_t i = first_nonfounder_index; i < event.size(); ++i) {
        stats[i] = static_cast<float>(event[i]);
    }
    //HOWTO: use std::copy with cast??    std::copy(event.begin()+first_nonfounder_index, event.end(), stats.begin() );
}


float MutationStats::get_mutation_prob() const {
    return mup;
}

bool MutationStats::get_has_single_mut() const {
    return has_single_mut;
}

const dng::GenotypeArray &MutationStats::inspect_posterior_at(int index) const {
    return posterior_probabilities[index];
}

const dng::GenotypeArray &MutationStats::inspect_genotype_at(int index) const {
    return genotype_likelihoods[index];
}


void MutationStats::record_basic_stats(hts::bcf::Variant &record){

    record.info("MUP", mup);
    record.info("LLD", lld);
    record.info("LLH", llh);
    record.info("MUX", mux); //OutputLevel::Complete
    record.info("MU1P", mu1p); //OutputLevel::single


};

void MutationStats::record_single_mutation_stats(hts::bcf::Variant &record){

    if(has_single_mut) {
        record.info("DNT", dnt);
        record.info("DNL", dnl);
        record.info("DNQ", dnq);
//        record.info("DNC", dnc);
        record.samples("MU1P", node_mu1p);
    }


};

void MutationStats::record_genotype_stats(hts::bcf::Variant &record){
    record.sample_genotypes(best_genotypes);
    record.samples("GQ", genotype_qualities);
    record.samples("GP", gp_scores);
    record.samples("GL", gl_scores);

}

void MutationStats::set_genotype_related_stats(const int (&acgt_to_refalt_allele)[5],
                                               const int (&refalt_to_acgt_allele)[5],
                                               const uint32_t n_alleles,
                                               const std::size_t ref_index,
                                               const std::size_t num_nodes,
                                               const std::size_t library_start
) {

// Construct numeric genotypes
    int numeric_genotype[10][2] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
                                   {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
    };
    for(int i = 0; i < 10; ++i) {
        int n1 = acgt_to_refalt_allele[dng::folded_diploid_nucleotides[i][0]];
        int n2 = acgt_to_refalt_allele[dng::folded_diploid_nucleotides[i][1]];
        if(n1 > n2) {
            numeric_genotype[i][0] = hts::bcf::encode_allele_unphased(n2);
            numeric_genotype[i][1] = hts::bcf::encode_allele_unphased(n1);
        } else {
            numeric_genotype[i][0] = hts::bcf::encode_allele_unphased(n1);
            numeric_genotype[i][1] = hts::bcf::encode_allele_unphased(n2);
        }
    }
    // Link VCF genotypes to our order
    int genotype_index[15];
    for(int i = 0, k = 0; i < n_alleles; ++i) {
        int n1 = refalt_to_acgt_allele[i];
        for(int j = 0; j <= i; ++j, ++k) {
            int n2 = refalt_to_acgt_allele[j];
            genotype_index[k] = (j == 0 && ref_index == 4) ?
                                -1 : dng::folded_diploid_genotypes_matrix[n1][n2];
        }
    }

    // Calculate sample genotypes
//    std::vector<int32_t> best_genotypes(2 * num_nodes);
//    std::vector<int32_t> genotype_qualities(num_nodes);
    int gt_count = n_alleles * (n_alleles + 1) / 2;
//    std::vector<float> gp_scores(gt_count*num_nodes );
    best_genotypes.resize(2 * num_nodes);
    genotype_qualities.resize(num_nodes);
    gp_scores.resize( gt_count * num_nodes);

    for(size_t i = 0, k = 0; i < num_nodes; ++i) {
        size_t pos;
        double d = posterior_probabilities[i].maxCoeff(&pos);
        best_genotypes[2 * i] = numeric_genotype[pos][0];
        best_genotypes[2 * i + 1] = numeric_genotype[pos][1];
        genotype_qualities[i] = dng::util::lphred<int32_t>(1.0 - d, 255);
        // If either of the alleles is missing set quality to 0
        if(hts::bcf::allele_is_missing({best_genotypes[2 * i]}) ||
                hts::bcf::allele_is_missing({best_genotypes[2 * i + 1]})) {
            genotype_qualities[i] = 0;
        }
        for(int j = 0; j < gt_count; ++j) {
            int n = genotype_index[j];
            gp_scores[k++] = (n == -1) ? 0.0 : posterior_probabilities[i][n];
        }
    }

    // Sample Likelihoods
//    std::vector<float> gl_scores.(gt_count *num_nodes, hts::bcf::float_missing);
    gl_scores.resize(gt_count *num_nodes, hts::bcf::float_missing);
    for(size_t i = library_start, k = library_start * gt_count; i < num_nodes;
        ++i) {
        for(int j = 0; j < gt_count; ++j) {
            int n = genotype_index[j];
            gl_scores[k++] = (n == -1) ? hts::bcf::float_missing :
                             genotype_likelihoods[i][n];
        }
    }


}