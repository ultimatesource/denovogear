/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
 * Copyright (c) 2015 Kael Dai
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *           Kael Dai <kdai1@asu.edu>
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

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <iosfwd>

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <sstream>
#include <string>

#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm/replace.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include <boost/algorithm/string.hpp>

#include <dng/task/call.h>
#include <dng/pedigree.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/likelihood.h>
#include <dng/seq.h>
#include <dng/utility.h>
#include <dng/hts/bcf.h>
#include <dng/hts/extra.h>
#include <dng/vcfpileup.h>
#include <dng/mutation.h>
#include <dng/stats.h>

#include <htslib/faidx.h>

#include "version.h"

using namespace dng;
using namespace dng::task;

// The main loop for dng-pileup application
// argument_type arg holds the processed command line arguments
int Pileup::operator()(Call::argument_type &arg) {
    using namespace std;

    // Open Reference
    io::Fasta reference{arg.fasta.c_str()};

    // quality thresholds
    int min_qual = arg.min_basequal;
    double min_prob = arg.min_prob;

    // Open input files
    vector<hts::bam::File> bamdata;
    for(auto && str : arg.input) {
        bamdata.emplace_back(str.c_str(), "r");
        if(bamdata.back().is_open()) {
            continue;
        }
        throw std::runtime_error("unable to open bam/sam file '" + str + "' for reading.");
    }

    // Add each genotype/sample column
    dng::ReadGroups rgs;
    rgs.ParseHeaderText(bamdata, arg.rgtag);

    // Pileup data
    std::vector<depth_t> read_depths(rgs.libraries().size());

    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    dng::BamPileup::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.qual.first[r.pos] < min_basequal
        || seq::base_index(r.aln.seq_at(r.pos)) >= 4);
    };

    // Read header from first file
    const bam_hdr_t *h = bamdata[0].header();
    int ref_target_id = -1, ref_sz = 0;
    const char* ref = nullptr;
    dng::BamPileup mpileup{rgs.groups(), arg.min_qlen};
    mpileup(bamdata, [ref_target_id, ref, ref_sz, h, &](const dng::BamPileup::data_type & data, uint64_t loc) {
        // Calculate target position and fetch sequence name
        int target_id = location_to_target(loc);
        int position = location_to_position(loc);

        if(target_id != ref_target_id && reference) {
            std::tie(ref, ref_sz) = reference.FetchSequence(h->target_name[target_id]);
            ref_target_id = target_id;
        }

        // Calculate reference base
        char ref_base = (ref && 0 <= position && position < ref_sz) ?
                        ref[position] : 'N';

        // reset all depth counters
        read_depths.assign(read_depths.size(), {});

        // pileup on read counts
        for(std::size_t u = 0; u < data.size(); ++u) {
            for(auto && r : data[u]) {
                if(filter_read(r)) {
                    continue;
                }
                std::size_t base = seq::base_index(r.aln.seq_at(r.pos));
                read_depths[rgs.library_from_id(u)].counts[ base ] += 1;
                // detect overflow of our signed number
                assert(read_depths[rgs.library_from_id(u)].counts[ base ] >= 0);
            }
        }

        // Determine what nucleotides show up and the order they will appear in the REF and ALT field
        // TODO: write tests that make sure REF="N" is properly handled
        //      (1) N should be included in AD only if REF="N"
        //      (2) N in AD should always be 0

        // Measure total depth and sort nucleotides in descending order
        typedef pair<int, int> key_t;
        key_t total_depths[4] = {{0, 0}, {1, 0}, {2, 0}, {3, 0}};
        int acgt_to_refalt_allele[5] = { -1, -1, -1, -1, -1}; // Maps allele to REF+ALT order
        int refalt_to_acgt_allele[5] = { -1, -1, -1, -1, -1}; // Maps REF+ALT order to A,C,G,T,N order

        int32_t dp_info = 0;
        std::vector<int32_t> dp_counts(num_nodes, hts::bcf::int32_missing);
        size_t dp_pos = library_start;
        for(auto && a : read_depths) {
            total_depths[0].second += a.counts[0];
            total_depths[1].second += a.counts[1];
            total_depths[2].second += a.counts[2];
            total_depths[3].second += a.counts[3];
            int32_t d = a.counts[0] + a.counts[1] + a.counts[2] + a.counts[3];
            dp_info += d;
            dp_counts[dp_pos++] = d;
        }
        sort(&total_depths[0], &total_depths[4], [](key_t a, key_t b) { return a.second > b.second; });

        // Construct a string representation of REF+ALT by ignoring nucleotides with no coverage
        string allele_order_str{seq::indexed_char(ref_index)};
        acgt_to_refalt_allele[ref_index] = 0;
        refalt_to_acgt_allele[0] = ref_index;
        int allele_count = 0; // Measures how many alleles in total_depths are non zero
        int refalt_count = 1; // Measures size of REF+ALT
        for(; allele_count < 4
                && total_depths[allele_count].second > 0; allele_count++) {
            if(total_depths[allele_count].first == ref_index) {
                continue;
            }
            allele_order_str += std::string(",") + seq::indexed_char(
                                    total_depths[allele_count].first);
            acgt_to_refalt_allele[total_depths[allele_count].first] = refalt_count;
            refalt_to_acgt_allele[refalt_count] = total_depths[allele_count].first;
            ++refalt_count;
        }
        // Update REF, ALT fields
        record.alleles(allele_order_str);

        // Construct numeric genotypes
        int numeric_genotype[10][2] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
        };
        for(int i = 0; i < 10; ++i) {
            int n1 = acgt_to_refalt_allele[folded_diploid_nucleotides[i][0]];
            int n2 = acgt_to_refalt_allele[folded_diploid_nucleotides[i][1]];
            if(n1 > n2) {
                numeric_genotype[i][0] = encode_allele_unphased(n2);
                numeric_genotype[i][1] = encode_allele_unphased(n1);
            } else {
                numeric_genotype[i][0] = encode_allele_unphased(n1);
                numeric_genotype[i][1] = encode_allele_unphased(n2);
            }
        }
        // Link VCF genotypes to our order
        int genotype_index[15];
        for(int i = 0, k = 0; i < refalt_count; ++i) {
            int n1 = refalt_to_acgt_allele[i];
            for(int j = 0; j <= i; ++j, ++k) {
                int n2 = refalt_to_acgt_allele[j];
                genotype_index[k] = (j == 0 && ref_index == 4) ?
                                    -1 : folded_diploid_genotypes_matrix[n1][n2];
            }
        }

        // Calculate sample genotypes
        vector<int32_t> best_genotypes(2 * num_nodes);
        vector<int32_t> genotype_qualities(num_nodes);
        int gt_count = refalt_count * (refalt_count + 1) / 2;
        vector<float> gp_scores(num_nodes * gt_count);

        for(size_t i = 0, k = 0; i < num_nodes; ++i) {
            size_t pos;
            double d = stats.posterior_probabilities[i].maxCoeff(&pos);
            best_genotypes[2 * i] = numeric_genotype[pos][0];
            best_genotypes[2 * i + 1] = numeric_genotype[pos][1];
            genotype_qualities[i] = lphred<int32_t>(1.0 - d, 255);
            // If either of the alleles is missing set quality to 0
            if(allele_is_missing({best_genotypes[2 * i]}) ||
                    allele_is_missing({best_genotypes[2 * i + 1]})) {
                genotype_qualities[i] = 0;
            }
            for(int j = 0; j < gt_count; ++j) {
                int n = genotype_index[j];
                gp_scores[k++] = (n == -1) ? 0.0 : stats.posterior_probabilities[i][n];
            }
        }

        // Sample Likelihoods
        vector<float> gl_scores(num_nodes * gt_count, hts::bcf::float_missing);
        for(size_t i = library_start, k = library_start * gt_count; i < num_nodes;
                ++i) {
            for(int j = 0; j < gt_count; ++j) {
                int n = genotype_index[j];
                gl_scores[k++] = (n == -1) ? hts::bcf::float_missing :
                                 stats.genotype_likelihoods[i][n];
            }
        }

        // Turn allele frequencies into AD format; order will need to match REF+ALT ordering of nucleotides
        vector<int32_t> ad_counts(num_nodes * refalt_count, hts::bcf::int32_missing);
        vector<int32_t> ad_info(refalt_count, 0);
        for(size_t u = 0; u < read_depths.size(); ++u) {
            size_t library_pos = (library_start + u) * refalt_count;
            for(size_t k = 0; k < refalt_count; ++k) {
                size_t index = refalt_to_acgt_allele[k];
                ad_counts[library_pos + k] = (index == 4) ? 0 : read_depths[u].counts[index];
                ad_info[k] += (index == 4) ? 0 : read_depths[u].counts[index];
            }
        }

        // Calculate ADR and ADF Tags
        vector<int32_t> adf_counts(num_nodes * refalt_count, hts::bcf::int32_missing);
        vector<int32_t> adr_counts(num_nodes * refalt_count, hts::bcf::int32_missing);
        vector<int32_t> adf_info(refalt_count, 0);
        vector<int32_t> adr_info(refalt_count, 0);

        vector<uint8_t> qual_ref, qual_alt;

        for(std::size_t u = library_start * refalt_count; u < num_nodes * refalt_count;
                ++u) {
            adf_counts[u] = 0;
            adr_counts[u] = 0;
        }
        for(std::size_t u = 0; u < data.size(); ++u) {
            std::size_t pos = (library_start + u) * refalt_count;
            for(auto && r : data[u]) {
                if(filter_read(r)) {
                    continue;
                }
                std::size_t library_pos = (library_start + rgs.library_from_id(
                                               u)) * refalt_count;
                std::size_t base = seq::base_index(r.aln.seq_at(r.pos));
                std::size_t base_refalt = acgt_to_refalt_allele[base];
                assert(base_refalt != -1);
                std::size_t base_pos = library_pos + base_refalt;
                // Forward Depths, avoiding branching
                adf_counts[base_pos]  += !r.aln.is_reversed();
                adf_info[base_refalt] += !r.aln.is_reversed();
                // Reverse Depths
                adr_counts[base_pos]  += r.aln.is_reversed();
                adr_info[base_refalt] += r.aln.is_reversed();
                // Mapping quality
                (base_refalt == 0 ? &qual_ref : &qual_alt)->push_back(r.aln.map_qual());
            }
        }

        // Fisher Exact Test for strand bias
        double fs_info;
        {
            int a11 = adf_info[0];
            int a21 = adr_info[0];
            int a12 = 0, a22 = 0;
            for(int k = 1; k < refalt_count; ++k) {
                a12 += adf_info[k];
                a22 += adr_info[k];
            }
            fs_info = dng::stats::fisher_exact_test(a11, a12, a21, a22);
        }
        double ta_info = dng::stats::ad_two_sample_test(qual_ref, qual_alt);

        record.info("MUP", stats.mup);
        record.info("LLD", stats.lld);
        record.info("LLH", stats.llh);
        record.info("MUX", stats.mux);
        record.info("MU1P", stats.mu1p);


        record.sample_genotypes(best_genotypes);
        record.samples("GQ", genotype_qualities);
        record.samples("GP", gp_scores);
        record.samples("GL", gl_scores);
        record.samples("DP", dp_counts);
        record.samples("AD", ad_counts);
        record.samples("ADF", adf_counts);
        record.samples("ADR", adr_counts);

        record.samples("MUP", stats.node_mup);

        if(stats.has_single_mut) {
            record.info("DNT", stats.dnt);
            record.info("DNL", stats.dnl);
            record.info("DNQ", stats.dnq);
            record.info("DNC", stats.dnc);

            record.samples("MU1P", stats.node_mu1p);
        }

        record.info("DP", dp_info);
        record.info("AD", ad_info);
        record.info("ADF", adf_info);
        record.info("ADR", adr_info);
        record.info("FS", static_cast<float>(phred(fs_info)));
        record.info("MQTa", static_cast<float>(ta_info));

        record.target(h->target_name[target_id]);
        record.position(position);
        vcfout.WriteRecord(record);
        record.Clear();
    });

    return EXIT_SUCCESS;
}

FindMutations::FindMutations(double min_prob, const Pedigree &pedigree,
                             params_t params) :
    pedigree_{pedigree}, min_prob_{min_prob},
    params_(params), genotype_likelihood_{params.params_a, params.params_b},
    work_(pedigree.CreateWorkspace()) {

    using namespace dng;

    // Use a parent-independent mutation model, which produces a
    // beta-binomial
    genotype_prior_[0] = population_prior(params_.theta, params_.nuc_freq, {params_.ref_weight, 0, 0, 0});
    genotype_prior_[1] = population_prior(params_.theta, params_.nuc_freq, {0, params_.ref_weight, 0, 0});
    genotype_prior_[2] = population_prior(params_.theta, params_.nuc_freq, {0, 0, params_.ref_weight, 0});
    genotype_prior_[3] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, params_.ref_weight});
    genotype_prior_[4] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, 0});

    // Calculate mutation expectation matrices
    full_transition_matrices_.assign(work_.num_nodes, {});
    nomut_transition_matrices_.assign(work_.num_nodes, {});
    posmut_transition_matrices_.assign(work_.num_nodes, {});
    onemut_transition_matrices_.assign(work_.num_nodes, {});
    mean_matrices_.assign(work_.num_nodes, {});

    for(size_t child = 0; child < work_.num_nodes; ++child) {
        auto trans = pedigree.transitions()[child];
        if(trans.type == Pedigree::TransitionType::Germline) {
            auto dad = f81::matrix(trans.length1, params_.nuc_freq);
            auto mom = f81::matrix(trans.length2, params_.nuc_freq);

            full_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom);
            nomut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 1);
            mean_matrices_[child] = meiosis_diploid_mean_matrix(dad, mom);
        } else if(trans.type == Pedigree::TransitionType::Somatic ||
                  trans.type == Pedigree::TransitionType::Library) {
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);

            full_transition_matrices_[child] = mitosis_diploid_matrix(orig);
            nomut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 1);
            mean_matrices_[child] = mitosis_diploid_mean_matrix(orig);
        } else {
            full_transition_matrices_[child] = {};
            nomut_transition_matrices_[child] = {};
            posmut_transition_matrices_[child] = {};
            onemut_transition_matrices_[child] = {};
            mean_matrices_[child] = {};
        }
    }

    //Calculate max_entropy based on having no data
    for(int ref_index = 0; ref_index < 5; ++ref_index) {
        work_.SetFounders(genotype_prior_[ref_index]);

        pedigree_.PeelForwards(work_, nomut_transition_matrices_);
        pedigree_.PeelBackwards(work_, nomut_transition_matrices_);
        event_.assign(work_.num_nodes, 0.0);
        double total = 0.0, entropy = 0.0;
        for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
            Eigen::ArrayXXd mat = (work_.super[i].matrix() *
                                   work_.lower[i].matrix().transpose()).array() *
                                  onemut_transition_matrices_[i].array();
            total += mat.sum();
            entropy += (mat.array() == 0.0).select(mat.array(),
                                                   mat.array() * mat.log()).sum();
        }
        // Calculate entropy of mutation location
        max_entropies_[ref_index] = (-entropy / total + log(total)) / M_LN2;
    }
}

// Returns true if a mutation was found and the record was modified
bool FindMutations::operator()(const std::vector<depth_t> &depths,
                               int ref_index, stats_t *stats) {
    using namespace std;
    using namespace hts::bcf;
    using dng::util::lphred;
    using dng::util::phred;

    assert(stats != nullptr);

    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::tie(work_.lower[work_.library_nodes.first + u], stemp) =
            genotype_likelihood_(depths[u], ref_index);
        scale += stemp;
    }

    // Set the prior probability of the founders given the reference
    work_.SetFounders(genotype_prior_[ref_index]);

    // Calculate log P(Data, nomut ; model)
    const double logdata_nomut = pedigree_.PeelForwards(work_,
                                 nomut_transition_matrices_);

    /**** Forward-Backwards with full-mutation ****/

    // Calculate log P(Data ; model)
    const double logdata = pedigree_.PeelForwards(work_, full_transition_matrices_);

    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    const double pmut = -std::expm1(logdata_nomut - logdata);

    // Skip this site if it does not meet lower probability threshold
    if(pmut < min_prob_) {
        return false;
    }

    // Copy genotype likelihoods
    stats->genotype_likelihoods.resize(work_.num_nodes);
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::size_t pos = work_.library_nodes.first + u;
        stats->genotype_likelihoods[pos] = work_.lower[pos].log() / M_LN10;
    }

    // Peel Backwards with full-mutation
    pedigree_.PeelBackwards(work_, full_transition_matrices_);

    size_t library_start = work_.library_nodes.first;

    stats->mup = pmut;
    stats->lld = (logdata + scale) / M_LN10;
    stats->llh = logdata / M_LN10;

    // Calculate statistics after Forward-Backwards
    stats->posterior_probabilities.resize(work_.num_nodes);
    for(std::size_t i = 0; i < work_.num_nodes; ++i) {
        stats->posterior_probabilities[i] = work_.upper[i] * work_.lower[i];
        stats->posterior_probabilities[i] /= stats->posterior_probabilities[i].sum();
    }
    double mux = 0.0;
    event_.assign(work_.num_nodes, 0.0);
    for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        mux += (work_.super[i] * (mean_matrices_[i] *
                                  work_.lower[i].matrix()).array()).sum();
        event_[i] = (work_.super[i] * (posmut_transition_matrices_[i] *
                                       work_.lower[i].matrix()).array()).sum();
        event_[i] = event_[i] / pmut;
    }
    stats->mux = mux;

    stats->node_mup.resize(work_.num_nodes, hts::bcf::float_missing);
    for(size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        stats->node_mup[i] = static_cast<float>(event_[i]);
    }

    /**** Forward-Backwards with no-mutation ****/

    // TODO: Better to use a separate workspace???
    pedigree_.PeelForwards(work_, nomut_transition_matrices_);
    pedigree_.PeelBackwards(work_, nomut_transition_matrices_);
    event_.assign(work_.num_nodes, 0.0);
    double total = 0.0, entropy = 0.0, max_coeff = -1.0;
    size_t dn_loc = 0, dn_col = 0, dn_row = 0;
    for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        Eigen::ArrayXXd mat = (work_.super[i].matrix() *
                               work_.lower[i].matrix().transpose()).array() *
                              onemut_transition_matrices_[i].array();
        std::size_t row, col;
        double mat_max = mat.maxCoeff(&row, &col);
        if(mat_max > max_coeff) {
            max_coeff = mat_max;
            dn_row  = row;
            dn_col = col;
            dn_loc = i;
        }
        event_[i] = mat.sum();
        entropy += (mat.array() == 0.0).select(mat.array(),
                                               mat.array() * mat.log()).sum();
        total += event_[i];
    }
    // Calculate P(only 1 mutation)
    const double pmut1 = total * (1.0 - pmut);
    stats->mu1p = pmut1;

    // Output statistics for single mutation only if it is likely
    if(pmut1 / pmut >= min_prob_) {
        stats->has_single_mut = true;
        for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
            event_[i] = event_[i] / total;
        }

        // Calculate entropy of mutation location
        entropy = (-entropy / total + log(total)) / M_LN2;
        entropy /= max_entropies_[ref_index];
        stats->dnc = std::round(100.0 * (1.0 - entropy));

        stats->dnq = lphred<int32_t>(1.0 - (max_coeff / total), 255);
        stats->dnl = pedigree_.labels()[dn_loc];
        if(pedigree_.transitions()[dn_loc].type == Pedigree::TransitionType::Germline) {
            stats->dnt = &meiotic_diploid_mutation_labels[dn_row][dn_col][0];
        } else {
            stats->dnt = &mitotic_diploid_mutation_labels[dn_row][dn_col][0];
        }

        stats->node_mu1p.resize(work_.num_nodes, hts::bcf::float_missing);
        for(size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
            stats->node_mu1p[i] = static_cast<float>(event_[i]);
        }
    } else {
        stats->has_single_mut = false;
    }
    return true;
}
