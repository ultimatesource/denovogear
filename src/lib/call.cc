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
#include <dng/utilities.h>
#include <dng/hts/bcf.h>
#include <dng/hts/extra.h>
#include <dng/vcfpileup.h>
#include <dng/mutation.h>
#include <dng/stats.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>

#include "version.h"

using namespace dng::task;
using namespace dng;

// Helper function that mimics boost::istream_range
template<class Elem, class Traits> inline
boost::iterator_range<std::istreambuf_iterator<Elem, Traits> >
istreambuf_range(std::basic_istream<Elem, Traits> &in) {
    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
               std::istreambuf_iterator<Elem, Traits>(in),
               std::istreambuf_iterator<Elem, Traits>());
}

// Helper function to determines if output should be bcf file, vcf file, or stdout. Also
// parses filename "bcf:<file>" --> "<file>"
std::pair<std::string, std::string> vcf_get_output_mode(
    Call::argument_type &arg) {
    using boost::algorithm::iequals;

    if(arg.output.empty() || arg.output == "-")
        return {"-", "w"};
    auto ret = hts::extra::extract_file_type(arg.output);
    if(iequals(ret.first, "bcf")) {
        return {ret.second, "wb"};
    } else if(iequals(ret.first, "vcf")) {
        return {ret.second, "w"};
    } else {
        throw std::runtime_error("Unknown file format '" + ret.second + "' for output '"
                                 + arg.output + "'.");
    }
    return {};
}

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
        str += std::string("--") + arg + '=' + dng::util::to_pretty(a) + ' ';
    }
    str.pop_back();
    return str;
}


template<typename VAL>
std::string vcf_command_line_text(const char *arg, VAL val) {
    return std::string("--") + arg + '=' + dng::util::to_pretty(val);
}

std::string vcf_command_line_text(const char *arg, std::string val) {
    return std::string("--") + arg + "=\'" + val + "\'";
}

// Helper function for writing the vcf header information
void vcf_add_header_text(hts::bcf::File &vcfout, Call::argument_type &arg) {
    using namespace std;
    string line{"##DeNovoGearCommandLine=<ID=dng-call,Version="
                PACKAGE_VERSION ","};
    line += vcf_timestamp();
    line += ",CommandLineOptions=\"";

#define XM(lname, sname, desc, type, def) \
	line += vcf_command_line_text(XS(lname),arg.XV(lname)) + ' ';
#	include <dng/task/call.xmh>
#undef XM
    for(auto && a : arg.input) {
        line += a + ' ';
    }

    line.pop_back();
    line += "\">";
    vcfout.AddHeaderMetadata(line);

    // Add the available tags for INFO, FILTER, and FORMAT fields
    vcfout.AddHeaderMetadata("##INFO=<ID=MUP,Number=1,Type=Float,Description=\"Probability of at least 1 de novo mutation\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=MU1P,Number=1,Type=Float,Description=\"Probability of exactly 1 de novo mutation\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=MUX,Number=1,Type=Float,Description=\"Expected number of de novo mutations\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=LLD,Number=1,Type=Float,Description=\"Log10-likelihood of observed data at the site\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=LLH,Number=1,Type=Float,Description=\"Scaled log10-likelihood of data at the site\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DNT,Number=1,Type=String,Description=\"De novo type\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DNL,Number=1,Type=String,Description=\"De novo location\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DNQ,Number=1,Type=Integer,Description=\"Phread-scaled de novo quality\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DNC,Number=1,Type=Integer,Description=\"De novo location certainty\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed (forward strand)\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed (reverse strand)\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=MQTa,Number=1,Type=Float,Description=\"Anderson-Darling Ta statistic for Alt vs. Ref read mapping qualities\">");

    vcfout.AddHeaderMetadata("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scaled genotype quality\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-likelihood of genotype based on read depths\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed (forward strand)\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed (reverse strand)\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=MUP,Number=1,Type=Float,Description=\"Probability of at least 1 de novo mutation in this node.\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=MU1P,Number=1,Type=Float,Description=\"Conditional probability that this node contains a de novo mutation given only 1 de novo mutation\">");
}

class FindMutations {
public:
    struct params_t {
        double theta;
        std::array<double, 4> nuc_freq;
        double ref_weight;

        dng::genotype::DirichletMultinomialMixture::params_t params_a;
        dng::genotype::DirichletMultinomialMixture::params_t params_b;
    };

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
};

// Build a list of all of the possible contigs to add to the vcf header
std::vector<std::pair<std::string, uint32_t>> parse_contigs(
const bam_hdr_t *hdr) {
    if(hdr == nullptr)
        return {};
    std::vector<std::pair<std::string, uint32_t>> contigs;
    uint32_t n_targets = hdr->n_targets;
    for(size_t a = 0; a < n_targets; a++) {
        if(hdr->target_name[a] == nullptr) {
            continue;
        }
        contigs.emplace_back(hdr->target_name[a], hdr->target_len[a]);
    }
    return contigs;
}


// VCF header lacks a function to get sequence lengths
// So we will extract the contig lines from the input header
std::vector<std::string> extract_contigs(const bcf_hdr_t *hdr) {
    if(hdr == nullptr)
        return {};
    // Read text of header
    int len;
    std::unique_ptr<char[], void(*)(void *)> str{bcf_hdr_fmt_text(hdr, 0, &len), free};
    if(!str)
        return {};
    std::vector<std::string> contigs;

    // parse ##contig lines
    const char *text = str.get();
    if(strncmp(text, "##contig=", 9) != 0) {
        text = strstr(text, "\n##contig=");
    } else {
        text = text - 1;
    }
    const char *end;
    for(; text != nullptr; text = strstr(end, "\n##contig=")) {
        for(end = text + 10; *end != '\n' && *end != '\0'; ++end)
            /*noop*/;
        if(*end != '\n') {
            return contigs;    // bad header, return what we have.
        }
        contigs.emplace_back(text + 1, end);
    }

    return contigs;
}

// The main loop for dng-call application
// argument_type arg holds the processed command line arguments
int Call::operator()(Call::argument_type &arg) {
    using namespace std;
    using namespace hts::bcf;
    using dng::util::lphred;
    using dng::util::phred;

    // Parse pedigree from file
    dng::io::Pedigree ped;
    if(!arg.ped.empty()) {
        ifstream ped_file(arg.ped);
        if(!ped_file.is_open()) {
            throw std::runtime_error(
                "unable to open pedigree file '" + arg.ped + "'.");
        }
        ped.Parse(istreambuf_range(ped_file));
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
    // TODO: read directly into freqs????  This will need a wrapper that provides an "insert" function.
    // TODO: include the size into the pattern, but this makes it harder to catch the second error.
    // TODO: turn all of this into a template function that returns array<double,4>?
    {
        auto f = util::parse_double_list(arg.nuc_freqs, ',', 4);
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
    double min_prob = arg.min_prob;

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

    // Begin writing VCF header
    auto out_file = vcf_get_output_mode(arg);
    hts::bcf::File vcfout(out_file.first.c_str(), out_file.second.c_str());
    vcf_add_header_text(vcfout, arg);

    if(cat == sequence_data) {
        // Wrap input in hts::bam::File
        for(auto && f : indata) {
            bamdata.emplace_back(std::move(f), arg.region.c_str(), arg.fasta.c_str(),
                                 arg.min_mapqual);
        }
        // Read header from first file
        const bam_hdr_t *h = bamdata[0].header();

        // Add contigs to header
        for(auto && contig : parse_contigs(h)) {
            vcfout.AddContig(contig.first.c_str(), contig.second);
        }

        // Add each genotype/sample column
        rgs.ParseHeaderText(bamdata);
    } else if(cat == variant_data) {
        bcfdata.emplace_back(std::move(indata[0]));
        // Read header from first file
        const bcf_hdr_t *h = bcfdata[0].header();

        // Add contigs to header
        for(auto && contig : extract_contigs(h)) {
            vcfout.AddHeaderMetadata(contig.c_str());
        }
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
        vcfout.AddHeaderMetadata(line.c_str());
    }

    FindMutations calculate { min_prob, pedigree,
        { arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1] } };

    // Pileup data
    std::vector<depth_t> read_depths(rgs.libraries().size());

    // Finish Header
    for(auto && str : pedigree.labels()) {
        vcfout.AddSample(str.c_str());
    }
    vcfout.WriteHeader();

    auto record = vcfout.InitVariant();
    const size_t num_nodes = pedigree.num_nodes();
    const size_t library_start = pedigree.library_nodes().first;

    FindMutations::stats_t stats;

    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    dng::BamPileup::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.qual.first[r.pos] < min_basequal
        || seq::base_index(r.aln.seq_at(r.pos)) >= 4);
    };

    // Treat sequence_data and variant data separately
    if(cat == sequence_data) {
        const bam_hdr_t *h = bamdata[0].header();
        dng::BamPileup mpileup{rgs.groups(), arg.min_qlen};
        mpileup(bamdata, [&](const dng::BamPileup::data_type & data, uint64_t loc) {

            // Calculate target position and fetch sequence name
            int target_id = location_to_target(loc);
            int position = location_to_position(loc);

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
            // TODO: handle overflow?  Down sample? Reservoir Sample?
            // The biggest read-depth we can handle per nucleotide is 2^16-1=65535
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
            if(!calculate(read_depths, ref_index, &stats)) {
                return;
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
            for(int a = 0; a < n_alleles; a++) {
                char base = *(rec->d.allele[a]);
                a2i.push_back(seq::char_index(base));
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

            if(!calculate(read_depths, seq::char_index(ref_base), &stats)) {
                return;
            }
            record.target(chrom);
            record.position(position);
            vcfout.WriteRecord(record);
            record.Clear();
        });
    } else {
        throw runtime_error("unsupported file category.");
    }

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
