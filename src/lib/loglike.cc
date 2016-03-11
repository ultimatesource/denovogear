/*
 * Copyright (c) 2016 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
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

#include <dng/task/loglike.h>
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
#include <dng/io/utility.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>

#include "version.h"

using namespace dng;

std::string vcf_timestamp() {
    using namespace std;
    using namespace std::chrono;
    std::string buffer(127, '\0');
    auto now = system_clock::now();
    auto now_t = system_clock::to_time_t(now);
    size_t sz = strftime(&buffer[0], 127, "Date=\"%FT%T%z\",Epoch=",
                         localtime(&now_t));
    buffer.resize(sz);
    auto epoch = std::chrono::duration_cast<std::chrono::milliseconds>(
                     now.time_since_epoch());
    buffer += to_string(epoch.count());
    return buffer;
}

template<typename V, typename A>
std::string vcf_command_line_text(const char *arg,
                                  const std::vector<V, A> &val) {
    std::string str;
    for(auto && a : val) {
        str += std::string("--") + arg + '=' + dng::utility::to_pretty(a) + ' ';
    }
    str.pop_back();
    return str;
}


template<typename VAL>
std::string vcf_command_line_text(const char *arg, VAL val) {
    return std::string("--") + arg + '=' + dng::utility::to_pretty(val);
}

std::string vcf_command_line_text(const char *arg, std::string val) {
    return std::string("--") + arg + "=\'" + val + "\'";
}

// Helper function for writing the vcf header information
void cout_add_header_text(task::LogLike::argument_type &arg) {
    using namespace std;
    string line{"##DeNovoGearCommandLine=<ID=dng-loglike,Version="
                PACKAGE_VERSION ","};
    line += vcf_timestamp();
    line += ",CommandLineOptions=\"";

#define XM(lname, sname, desc, type, def) \
    line += vcf_command_line_text(XS(lname),arg.XV(lname)) + ' ';
#   include <dng/task/loglike.xmh>
#undef XM
    for(auto && a : arg.input) {
        line += a + ' ';
    }

    line.pop_back();
    line += "\">";

    std::cout << line << "\n";
 }

class LogProbability {
public:
    struct params_t {
        double theta;
        std::array<double, 4> nuc_freq;
        double ref_weight;

        dng::genotype::DirichletMultinomialMixture::params_t params_a;
        dng::genotype::DirichletMultinomialMixture::params_t params_b;
    };

    struct stats_t {
        double log_data;
        double log_scale;        
    };

    LogProbability(const Pedigree &pedigree, params_t params);

    stats_t operator()(const std::vector<depth_t> &depths, int ref_index);

protected:
    const dng::Pedigree &pedigree_;

    params_t params_;

    dng::peel::workspace_t work_;


    dng::TransitionVector full_transition_matrices_;

    // Model genotype likelihoods as a mixture of two dirichlet multinomials
    // TODO: control these with parameters
    dng::genotype::DirichletMultinomialMixture genotype_likelihood_;

    dng::GenotypeArray genotype_prior_[5]; // Holds P(G | theta)
};

// The main loop for dng-loglike application
// argument_type arg holds the processed command line arguments
int task::LogLike::operator()(task::LogLike::argument_type &arg) {
    using namespace std;
    using namespace hts::bcf;

    // Parse pedigree from file
    dng::io::Pedigree ped;
    if(!arg.ped.empty()) {
        ifstream ped_file(arg.ped);
        if(!ped_file.is_open()) {
            throw std::runtime_error(
                "unable to open pedigree file '" + arg.ped + "'.");
        }
        ped.Parse(io::istreambuf_range(ped_file));
    } else {
        throw std::runtime_error("pedigree file was not specified.");
    }

    // Open Reference
    unique_ptr<char[], void(*)(void *)> ref{nullptr, free};
    int ref_sz = 0, ref_target_id = -1;
    unique_ptr<faidx_t, void(*)(faidx_t *)> fai{nullptr, fai_destroy};
    if(!arg.fasta.empty()) {
        fai.reset(fai_load(arg.fasta.c_str()));
        if(!fai)
            throw std::runtime_error("unable to open faidx-indexed reference file '"
                                     + arg.fasta + "'.");
    }

    // Parse Nucleotide Frequencies
    std::array<double, 4> freqs;
    {
        auto f = utility::parse_double_list(arg.nuc_freqs, ',', 4);
        if(!f.second) {
            throw std::runtime_error("Unable to parse nuc-freq option. "
                                     "It must be a comma separated list of floating-point numbers.");
        }
        if(f.first.size() != 4) {
            throw std::runtime_error("Wrong number of values passed to nuc-freq. "
                                     "Expected 4; found " + std::to_string(f.first.size()) + ".");
        }
        std::copy(f.first.begin(), f.first.end(), &freqs[0]);
    }

    // quality thresholds
    int min_qual = arg.min_basequal;

    // Open input files
    dng::ReadGroups rgs;
    vector<hts::File> indata;
    vector<hts::bam::File> bamdata;
    vector<hts::bcf::File> bcfdata;
    for(auto && str : arg.input) {
        indata.emplace_back(str.c_str(), "r");
        if(indata.back().is_open()) {
            continue;
        }
        throw std::runtime_error("unable to open input file '" + str + "'.");
    }

    // Check to see if all inputs are of the same type
    const htsFormatCategory cat = indata[0].format().category;
    for(auto && f : indata) {
        if(f.format().category == cat) {
            continue;
        }
        throw std::runtime_error("mixing sequence data and variant data as input is not supported.");
    }

    cout_add_header_text(arg);

    // replace arg.region with the contents of a file if needed
    io::at_slurp(arg.region);

    if(cat == sequence_data) {
        // Wrap input in hts::bam::File
        for(auto && f : indata) {
       	    bamdata.emplace_back(std::move(f), arg.fasta.c_str(),
				 arg.min_mapqual, arg.header.c_str());
        }
        if(!arg.region.empty()) {
            for(auto && f : bamdata) {
                auto r = regions::bam_parse(arg.region, f);
                f.regions(std::move(r));
            }
        }        

        // Add each genotype/sample column
        rgs.ParseHeaderText(bamdata, arg.rgtag);
    } else if(cat == variant_data) {
        bcfdata.emplace_back(std::move(indata[0]));

        // Add each genotype/sample column
        rgs.ParseSamples(bcfdata[0]);
    } else {
        throw runtime_error("unsupported file category.");
    }

    // Construct peeling algorithm from parameters and pedigree information
    dng::Pedigree pedigree;
    if(!pedigree.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop pedigree.");
    }
    if(arg.gamma.size() < 2) {
        throw std::runtime_error("Unable to construct genotype-likelihood model; "
                                 "Gamma needs to be specified at least twice to change model from default.");
    }

    for(auto && line : pedigree.BCFHeaderLines()) {
        std:cout << line << "\n";
    }

    LogProbability calculate (pedigree,
        { arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1] } );

    // Pileup data
    std::vector<depth_t> read_depths(rgs.libraries().size());

    const size_t num_nodes = pedigree.num_nodes();
    const size_t library_start = pedigree.library_nodes().first;

    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    dng::BamPileup::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.qual.first[r.pos] < min_basequal
        || seq::base_index(r.aln.seq_at(r.pos)) >= 4);
    };

    // Treat sequence_data and variant data separately
    dng::stats::ExactSum sum_data;
    dng::stats::ExactSum sum_scale;
    if(cat == sequence_data) {
        const bam_hdr_t *h = bamdata[0].header();
        dng::BamPileup mpileup{rgs.groups(), arg.min_qlen};
        mpileup(bamdata, [&](const dng::BamPileup::data_type & data, location_t loc) {

            // Calculate target position and fetch sequence name
            int target_id = utility::location_to_target(loc);
            int position = utility::location_to_position(loc);

            if(target_id != ref_target_id && fai) {
                ref.reset(faidx_fetch_seq(fai.get(), h->target_name[target_id],
                                          0, 0x7fffffff, &ref_sz));
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
                    assert(read_depths[rgs.library_from_id(u)].counts[ base ] < 65535);

                    read_depths[rgs.library_from_id(u)].counts[ base ] += 1;
                }
            }
            size_t ref_index = seq::char_index(ref_base);
            auto loglike = calculate(read_depths, ref_index);
            sum_data += loglike.log_data;
            sum_scale += loglike.log_scale;
        });
    } else if(cat == variant_data) {
        const char *fname = bcfdata[0].name();
        dng::pileup::vcf::VCFPileup vcfpileup{rgs.libraries()};
        vcfpileup(fname, [&](bcf_hdr_t *hdr, bcf1_t *rec) {
            // Won't be able to access ref->d unless we unpack the record first
            bcf_unpack(rec, BCF_UN_STR);

            // get chrom, position, ref from their fields
            const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
            int32_t position = rec->pos;
            uint32_t n_alleles = rec->n_allele;
            uint32_t n_samples = bcf_hdr_nsamples(hdr);
            const char ref_base = *(rec->d.allele[0]);

            // Read all the Allele Depths for every sample into ad array
            int *ad = NULL;
            int n_ad = 0;
            int n_ad_array = 0;
            n_ad = bcf_get_format_int32(hdr, rec, "AD", &ad, &n_ad_array);

            // Create a map between the order of vcf alleles (REF+ALT) and their correct index in read_depths.counts[]
            vector<size_t> a2i;
            string allele_order_str;
            int acgt_to_refalt_allele[5] = { -1, -1, -1, -1, -1}; // Maps allele to REF+ALT order
            int refalt_to_acgt_allele[5] = { -1, -1, -1, -1, -1}; // Maps REF+ALT order to A,C,G,T,N order
            for(int a = 0; a < n_alleles; a++) {
                char base = *(rec->d.allele[a]);
                size_t base_indx = seq::char_index(base);
                a2i.push_back(seq::char_index(base));
                acgt_to_refalt_allele[base_indx] = a;
                refalt_to_acgt_allele[a] = base_indx;

                if(a != 0)
            		allele_order_str += ",";
            	allele_order_str += rec->d.allele[a];
            }

            // Build the read_depths
            read_depths.assign(n_samples, {});
            for(size_t sample_ndx = 0; sample_ndx < n_samples; sample_ndx++) {
                size_t pos = rgs.library_from_index(sample_ndx);
                if(pos == -1) {
                    continue;
                }
                for(size_t allele_ndx = 0; allele_ndx < n_alleles; allele_ndx++) {
                    int32_t depth = ad[n_alleles * sample_ndx + allele_ndx];
                    size_t base = a2i[allele_ndx];
                    if(!(base < 4)) {
                        continue;
                    }
                    read_depths[pos].counts[base] = depth;
                }
            }

            size_t ref_index = seq::char_index(ref_base);
            auto loglike = calculate(read_depths, ref_index);
            sum_data += loglike.log_data;
            sum_scale += loglike.log_scale;
        });
    } else {
        throw runtime_error("unsupported file category.");
    }
    std::cout << "log_hidden\tlog_observed\n";
    std::cout << setprecision(16) << sum_data.result() << "\t" << sum_scale.result() << "\n";

    return EXIT_SUCCESS;
}

LogProbability::LogProbability(const Pedigree &pedigree,
                             params_t params) :
    pedigree_(pedigree),
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
 
    for(size_t child = 0; child < work_.num_nodes; ++child) {
        auto trans = pedigree.transitions()[child];
        if(trans.type == Pedigree::TransitionType::Germline) {
            auto dad = f81::matrix(trans.length1, params_.nuc_freq);
            auto mom = f81::matrix(trans.length2, params_.nuc_freq);

            full_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom);
        } else if(trans.type == Pedigree::TransitionType::Somatic ||
                  trans.type == Pedigree::TransitionType::Library) {
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);

            full_transition_matrices_[child] = mitosis_diploid_matrix(orig);
        } else {
            full_transition_matrices_[child] = {};
        }
    }
}

// Returns true if a mutation was found and the record was modified
LogProbability::stats_t LogProbability::operator()(const std::vector<depth_t> &depths,
                               int ref_index) {
    using namespace std;

    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::tie(work_.lower[work_.library_nodes.first + u], stemp) =
            genotype_likelihood_(depths[u], ref_index);
        scale += stemp;
    }

    // Set the prior probability of the founders given the reference
    work_.SetFounders(genotype_prior_[ref_index]);

    // Calculate log P(Data ; model)
    double logdata = pedigree_.PeelForwards(work_, full_transition_matrices_);

    return {logdata/M_LN10, scale/M_LN10};
}
