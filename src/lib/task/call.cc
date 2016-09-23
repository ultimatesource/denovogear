/*
 * Copyright (c) 2014-2016 Reed A. Cartwright
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

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <sstream>
#include <string>

#include <boost/range/algorithm/replace.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include <boost/algorithm/string.hpp>

#include <dng/task/call.h>
#include <dng/relationship_graph.h>
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
#include <dng/find_mutations.h>

#include "htslib/synced_bcf_reader.h"

#include <htslib/faidx.h>
#include <htslib/khash.h>

#include "version.h"

using namespace dng;

// Helper function to determines if output should be bcf file, vcf file, or stdout. Also
// parses filename "bcf:<file>" --> "<file>"
std::pair<std::string, std::string> vcf_get_output_mode(
    task::Call::argument_type &arg) {
    using boost::algorithm::iequals;

    if(arg.output.empty() || arg.output == "-")
        return {"-", "w"};
    auto ret = utility::extract_file_type(arg.output);
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
void vcf_add_header_text(hts::bcf::File &vcfout, task::Call::argument_type &arg) {
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
    vcfout.AddHeaderMetadata("##INFO=<ID=LLS,Number=1,Type=Float,Description=\"Scaled log10-likelihood of observed data at the site\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DNT,Number=1,Type=String,Description=\"De novo type\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DNL,Number=1,Type=String,Description=\"De novo location\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DNQ,Number=1,Type=Integer,Description=\"Phread-scaled de novo quality\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DNC,Number=1,Type=Integer,Description=\"De novo location certainty\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed (forward strand)\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed (reverse strand)\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=MQTa,Number=1,Type=Float,Description=\"Anderson-Darling Ta statistic for Alt vs. Ref read mapping qualities\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=RPTa,Number=1,Type=Float,Description=\"Anderson-Darling Ta statistic for Alt vs. Ref read positions\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=BQTa,Number=1,Type=Float,Description=\"Anderson-Darling Ta statistic for Alt vs. Ref base-call qualities\">");

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



// The main loop for dng-call application
// argument_type arg holds the processed command line arguments
int task::Call::operator()(Call::argument_type &arg) {
    using namespace std;
    using namespace hts::bcf;
    using dng::utility::lphred;
    using dng::utility::phred;
    using dng::utility::location_to_contig;
    using dng::utility::location_to_position;

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
    // TODO: read directly into freqs????  This will need a wrapper that provides an "insert" function.
    // TODO: include the size into the pattern, but this makes it harder to catch the second error.
    // TODO: turn all of this into a template function that returns array<double,4>?
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
    // Scale frequencies and makes sure they sum to 1.
    for(int i=0;i<freqs.size();++i) {
        if(freqs[i] < 0.0) {
            throw std::runtime_error("Nucleotide frequencies must be non-negative. "
                "Parsing of '" + arg.nuc_freqs + "' generated a frequency of " + std::to_string(freqs[i])
                + " in position " + std::to_string(i) + "."
                );
        }
    }
    double freqs_sum = freqs[0] + freqs[1] + freqs[2] + freqs[3];
    for(int i=0;i<freqs.size();++i) {
        freqs[i] = freqs[i]/freqs_sum;
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

    // replace arg.region with the contents of a file if needed
    io::at_slurp(arg.region);

    bcf_srs_t *rec_reader;

    InheritanceModel inheritance_model;
    inheritance_model.parse_model(arg.model);

    if(cat == sequence_data) {
        // Wrap input in hts::bam::File
        for(auto && f : indata) {
       	    bamdata.emplace_back(std::move(f), arg.fasta.c_str(),
				 arg.min_mapqual, arg.header.c_str());
        }

        if(!arg.region.empty()) {
			if(arg.region.find(".bed")!=std::string::npos){
				for(auto && f : bamdata) {
					auto r = regions::bam_parse_bed(arg.region, f);
					f.regions(std::move(r));
				}
			}else{
				for(auto && f : bamdata) {
					auto r = regions::bam_parse_region(arg.region, f);
					f.regions(std::move(r));
				}
			}
        }

        // Read header from first file
        const bam_hdr_t *h = bamdata[0].header();

        // Add contigs to header
        for(auto && contig : hts::extra::parse_contigs(h)) {
            vcfout.AddContig(contig.first.c_str(), contig.second);
        }

        // Add each genotype/sample column
        rgs.ParseHeaderText(bamdata, arg.rgtag);
    } else if(cat == variant_data) {
        bcfdata.emplace_back(std::move(indata[0]));
        rec_reader = bcf_sr_init(); // used for iterating each rec in BCF/VCF

		// Open region if specified
		if(!arg.region.empty()) {
			int is_file = (arg.region.find("bed") != std::string::npos)? 1 : 0;
			int ret = bcf_sr_set_regions(rec_reader, arg.region.c_str(), is_file);
			if(ret == -1) {
				throw std::runtime_error("no records in the query region " + arg.region);
			}
		}

		// Initialize the record reader to iterate through the BCF/VCF input
		int ret = bcf_sr_add_reader(rec_reader, bcfdata[0].name());
		if(ret == 0) {
			int errnum = rec_reader->errnum;
			switch(errnum) {
			case not_bgzf:
				throw std::runtime_error("Input file type does not allow for region searchs. Exiting!");
				break;
			case idx_load_failed:
				throw std::runtime_error("Unable to load query region, no index. Exiting!");
				break;
			case file_type_error:
				throw std::runtime_error("Could not load filetype. Exiting!");
				break;
			default:
				throw std::runtime_error("Could not load input sequence file into htslib. Exiting!");
			};
		}

        // Read header from first file
        const bcf_hdr_t *h = bcf_sr_get_header(rec_reader, 0);

        // Add contigs to header
        for(auto && contig : hts::extra::extract_contigs(h)) {
            vcfout.AddHeaderMetadata(contig.c_str());
        }
        // Add each genotype/sample column
        rgs.ParseSamples(bcfdata[0]);
    } else {
        throw runtime_error("unsupported file category.");
    }

    // Construct peeling algorithm from parameters and pedigree information
    dng::RelationshipGraph relationship_graph;
    if (!relationship_graph.Construct(ped, rgs,
                                      inheritance_model.GetInheritancePattern(),
                                      arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop relationship_graph.");
    }
    if(arg.gamma.size() < 2) {
        throw std::runtime_error("Unable to construct genotype-likelihood model; "
                                 "Gamma needs to be specified at least twice to change model from default.");
    }

    for(auto && line : relationship_graph.BCFHeaderLines()) {
        vcfout.AddHeaderMetadata(line.c_str());
    }

    FindMutations calculate ( min_prob, relationship_graph,
        { arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1] } );

    // Pileup data
    std::vector<depth_t> read_depths(rgs.libraries().size());

    // Finish Header
    for(auto && str : relationship_graph.labels()) {
        vcfout.AddSample(str.c_str());
    }
    vcfout.WriteHeader();

    auto record = vcfout.InitVariant();
    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;

    FindMutations::stats_t stats;

    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    dng::BamPileup::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.base_qual() < min_basequal
        || seq::base_index(r.base()) >= 4);
    };

    // Treat sequence_data and variant data separately
    if(cat == sequence_data) {
        const bam_hdr_t *h = bamdata[0].header();
        dng::BamPileup mpileup{rgs.groups(), arg.min_qlen};
        mpileup(bamdata, [&](const dng::BamPileup::data_type & data, utility::location_t loc) {

            // Calculate target position and fetch sequence name
            int target_id = location_to_contig(loc);
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
            // The biggest read-depth we can handle per nucleotide is 2^32-1
            for(std::size_t u = 0; u < data.size(); ++u) {
                for(auto && r : data[u]) {
                    if(filter_read(r)) {
                        continue;
                    }
                    std::size_t base = seq::base_index(r.base());
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
            // calculate the log-likelihood of the null hypothesis that all reads come from binomial
            double log_null = 0.0;
            if(dp_info > 0.0) {
                log_null = -dp_info*log10(dp_info);
                for(int i=0;i<4;++i) {
                    if(total_depths[i].second > 0) {
                        log_null += total_depths[i].second*log10(total_depths[i].second);
                    }
                }
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
                    genotype_index[k] = (n1 == 4 || n2 == 4) ?
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

            vector<int> qual_ref, qual_alt, pos_ref, pos_alt, base_ref, base_alt;
            qual_ref.reserve(data.size());
            qual_alt.reserve(data.size());
            pos_ref.reserve(data.size());
            pos_alt.reserve(data.size());
            base_ref.reserve(data.size());
            base_alt.reserve(data.size());
            double rms_mq = 0.0;

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
                    std::size_t base = seq::base_index(r.base());
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
                    rms_mq += r.aln.map_qual()*r.aln.map_qual();
                    (base_refalt == 0 ? &qual_ref : &qual_alt)->push_back(r.aln.map_qual());
                    // Positions
                    (base_refalt == 0 ? &pos_ref : &pos_alt)->push_back(r.pos);
                    // Base Calls
                    (base_refalt == 0 ? &base_ref : &base_alt)->push_back(r.base_qual());
                }
            }
            rms_mq = sqrt(rms_mq/(qual_ref.size()+qual_alt.size()));

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
            double mq_info = dng::stats::ad_two_sample_test(qual_ref, qual_alt);
            double rp_info = dng::stats::ad_two_sample_test(pos_ref, pos_alt);
            double bq_info = dng::stats::ad_two_sample_test(base_ref, base_alt);

            record.info("MUP", stats.mup);
            record.info("LLD", stats.lld);
            record.info("LLS", static_cast<float>(stats.lld-log_null));
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
            record.info("MQ", static_cast<float>(rms_mq));
            record.info("FS", static_cast<float>(phred(fs_info)));
            record.info("MQTa", static_cast<float>(mq_info));
            record.info("RPTa", static_cast<float>(rp_info));
            record.info("BQTa", static_cast<float>(bq_info));

            record.target(h->target_name[target_id]);
            record.position(position);
            vcfout.WriteRecord(record);
            record.Clear();
        });
    } else if(cat == variant_data) {
        const char *fname = bcfdata[0].name();
        dng::pileup::vcf::VCFPileup vcfpileup{rgs.libraries()};
        vcfpileup(fname, rec_reader, [&](bcf_hdr_t *hdr, bcf1_t *rec){
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
            if(!calculate(read_depths, ref_index, &stats)) {
                return;
            }

            // reformatted AD fields for output. The first few fields are the GL and SM fields and "AD" is missing. The
            // remaining fields are just copied from the input file
            std::vector<int32_t> ad_counts(library_start*n_alleles, hts::bcf::int32_missing);
            ad_counts.insert(ad_counts.end(), ad, ad + n_samples*n_alleles);

            // sum up all the counts
            vector<int32_t> ad_info(n_alleles, 0);
            for(size_t a = 0; a < n_ad; a++)
            	ad_info[a%n_alleles] += ad[a];

            // sum up the depths for each sample and over all
            std::vector<int32_t> dp_counts(num_nodes, hts::bcf::int32_missing);
            int32_t dp_info = 0;
            int32_t total_depths[4] = {0,0,0,0};
            for(int sample = 0; sample < n_samples; sample++) {
            	int32_t count = 0;
            	for(int allele = 0; allele < n_alleles; allele++) {
            		count += ad[sample*n_alleles + allele];
                    total_depths[allele] += ad[sample*n_alleles + allele];
            	}
            	dp_counts[sample+library_start] = count;
            	dp_info += count;
            }
            // calculate the log-likelihood of the null hypothesis that all reads come from binomial
            double log_null = 0.0;
            if(dp_info > 0.0) {
                log_null = -dp_info*log10(dp_info);
                for(int allele=0;allele<n_alleles;++allele) {
                    if(total_depths[allele] > 0) {
                        log_null += total_depths[allele]*log10(total_depths[allele]);
                    }
                }
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
            for(int i = 0, k = 0; i < n_alleles; ++i) {
                int n1 = refalt_to_acgt_allele[i];
                for(int j = 0; j <= i; ++j, ++k) {
                    int n2 = refalt_to_acgt_allele[j];
                    genotype_index[k] = (n1 == 4 || n2 == 4) ?
                                        -1 : folded_diploid_genotypes_matrix[n1][n2];
                }
            }

            // Calculate sample genotypes
            vector<int32_t> best_genotypes(2 * num_nodes);
            vector<int32_t> genotype_qualities(num_nodes);
            int gt_count = n_alleles * (n_alleles + 1) / 2;
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


            record.info("MUP", stats.mup);
            record.info("LLD", stats.lld);
            record.info("LLS", static_cast<float>(stats.lld-log_null));
            record.info("MUX", stats.mux);
            record.info("MU1P", stats.mu1p);

            record.sample_genotypes(best_genotypes);
            record.samples("GQ", genotype_qualities);
            record.samples("GP", gp_scores);
            record.samples("GL", gl_scores);
            record.samples("DP", dp_counts);
            record.samples("AD", ad_counts);

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
