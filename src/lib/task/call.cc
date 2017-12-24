/*
 * Copyright (c) 2014-2016 Reed A. Cartwright
 * Copyright (c) 2015-2016 Kael Dai
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *           Kael Dai <kdai1@asu.edu>
 *           Steven H. Wu <stevenwu@asu.edu>
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
#include <boost/range/algorithm/fill.hpp>
#include <boost/range/algorithm/replace_if.hpp>

#include <boost/algorithm/string.hpp>

#include <dng/task/call.h>
#include <dng/relationship_graph.h>
#include <dng/fileio.h>
#include <dng/seq.h>
#include <dng/utility.h>
#include <dng/mutation.h>
#include <dng/stats.h>
#include <dng/io/utility.h>
#include <dng/io/fasta.h>
#include <dng/call_mutations.h>

#include <dng/io/bam.h>
#include <dng/io/bcf.h>
#include <dng/io/ad.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>

#include "version.h"

using namespace dng;

using utility::make_array;

int process_bam(task::Call::argument_type &arg);
int process_bcf(task::Call::argument_type &arg);
int process_ad(task::Call::argument_type &arg);

void add_stats_to_output(const CallMutations::stats_t& call_stats, const pileup::stats_t& depth_stats,
    bool has_single_mut,
    const RelationshipGraph &graph,
    const peel::workspace_t &work,
    hts::bcf::Variant *record);

// Helper function to determines if output should be bcf file, vcf file, or stdout. Also
// parses filename "bcf:<file>" --> "<file>"
std::pair<std::string, std::string> vcf_get_output_mode(
    const task::Call::argument_type &arg) {
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

// Helper function for writing the vcf header information
void vcf_add_header_text(hts::bcf::File &vcfout, const task::Call::argument_type &arg) {
    using namespace std;
    string line{"##DeNovoGearCommandLine=<ID=dng-call,Version="
                PACKAGE_VERSION ","};
    line += utility::vcf_timestamp();
    line += ",CommandLineOptions=\"";

#define XM(lname, sname, desc, type, def) \
	line +=  utility::vcf_command_line_text(XS(lname),arg.XV(lname)) + ' ';
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
    vcfout.AddHeaderMetadata("##FORMAT=<ID=MUP,Number=1,Type=Float,Description=\"Conditional probability that this node contains at least 1 de novo mutation given at least one mutation at this site\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=MU1P,Number=1,Type=Float,Description=\"Conditional probability that this node contains a de novo mutation given only 1 de novo mutation at this site\">");
}

template<typename A, typename M, typename R>
hts::bcf::File open_vcf_output(const A& arg, const M& mpileup, const R& relationship_graph) {
   // Begin writing VCF header
    auto out_file = vcf_get_output_mode(arg);
    hts::bcf::File vcfout(out_file.first.c_str(), out_file.second.c_str());
    vcf_add_header_text(vcfout, arg);

    for(auto && contig : mpileup.contigs()) {
        vcfout.AddContig(contig.name.c_str(), contig.length); // Add contigs to header
    }
    for(auto && line : relationship_graph.BCFHeaderLines()) {
        vcfout.AddHeaderMetadata(line.c_str());  // Add pedigree info
    }

    // Finish Header
    for(auto && str : relationship_graph.labels()) {
        vcfout.AddSample(str.c_str()); // Add genotype columns
    }
    vcfout.WriteHeader();

    return vcfout;
}

using namespace hts::bcf;
using dng::utility::lphred;
using dng::utility::phred;
using dng::utility::location_to_contig;
using dng::utility::location_to_position;
using dng::utility::FileCat;
using namespace task;

// The main loop for dng-call application
// argument_type arg holds the processed command line arguments
int task::Call::operator()(Call::argument_type &arg) {
    // Determine the type of input files
    auto it = arg.input.begin();
    FileCat mode = utility::input_category(*it, FileCat::Sequence|FileCat::Pileup|FileCat::Variant, FileCat::Unknown);
    for(++it; it != arg.input.end(); ++it) {
    	// Make sure different types of input files aren't mixed together
        if(utility::input_category(*it, FileCat::Sequence|FileCat::Pileup|FileCat::Variant, FileCat::Sequence) != mode) {
            throw std::invalid_argument("Mixing pileup, sequencing, and variant file types is not supported.");
        }
    }

    if(mode == utility::FileCat::Sequence) {
    	// sam, bam, cram
    	return process_bam(arg);
    } else if(mode == utility::FileCat::Variant) {
    	// vcf, bcf
    	return process_bcf(arg);
    } else if(mode == utility::FileCat::Pileup) {
    	// tad, ad
    	return process_ad(arg);
    } else {
    	throw std::invalid_argument("Unknown input data file type.");
    }
    return EXIT_FAILURE;
}

// Processes bam, sam, and cram files.
int process_bam(task::Call::argument_type &arg) {
    // Open Reference
    if(arg.fasta.empty()){
    	throw std::invalid_argument("Path to reference file must be specified with --fasta when processing bam/sam/cram files.");
    }
    io::Fasta reference{arg.fasta.c_str()};

    // Open input files
    auto mpileup = io::BamPileup::open_and_setup(arg);

    auto relationship_graph = create_relationship_graph(arg, &mpileup);

    // Open Output
    auto vcfout = open_vcf_output(arg, mpileup, relationship_graph);

    // Record for each output
    auto record = vcfout.InitVariant();

    // Construct Calling Object
    CallMutations do_call(arg.min_prob, relationship_graph, get_model_parameters(arg));

    // Pileup data
    dng::pileup::allele_depths_t read_depths(make_array(mpileup.num_libraries(),5u));

    // Calculated stats
    CallMutations::stats_t stats;
  
    // Parameters used by site calculation function
    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;
    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    decltype(mpileup)::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.base_qual() < min_basequal
        || seq::base_index(r.base()) >= 4);
    };

    auto h = mpileup.header(); // TODO: Fixme.

    const double min_prob = arg.min_prob;
    mpileup([&](const decltype(mpileup)::data_type & data, utility::location_t loc) {
        // Calculate target position and fetch sequence name
        int contig = utility::location_to_contig(loc);
        int position = utility::location_to_position(loc);

        // Calculate reference base
        assert(0 <= contig && contig < h->n_targets);
        char ref_base = reference.FetchBase(h->target_name[contig],position);
        size_t ref_index = seq::char_index(ref_base);

        // reset all depth counters
        std::fill_n(read_depths.data(), read_depths.num_elements(), 0);

        // pileup on read counts
        for(std::size_t u = 0; u < data.size(); ++u) {
            for(auto && r : data[u]) {
                if(filter_read(r)) {
                    continue;
                }
                std::size_t base = seq::base_index(r.base());
                assert(read_depths[u][(ref_index+base) % 5] < 65535);

                read_depths[u][(ref_index+base) % 5] += 1;
            }
        }
        int num_alts = (ref_index < 4) ? 3 : 4;

		if(!do_call(read_depths, num_alts, ref_index < 4, &stats)) {
			return;
		}
        const bool has_single_mut = ((stats.mu1p / stats.mup) >= min_prob);

		// Measure total depth and sort nucleotides in descending order
        pileup::stats_t depth_stats;
        pileup::calculate_stats(read_depths, &depth_stats);

        add_stats_to_output(stats, depth_stats, has_single_mut, relationship_graph, do_call.work(), &record);

        // Information about the site, based on depths
        const int color = 40;//depth_stats.color;
        const auto & type_info_gt = dng::pileup::AlleleDepths::type_info_gt_table[color];
        const auto & type_info = dng::pileup::AlleleDepths::type_info_table[color];
        const int refalt_count = type_info.width + (type_info.reference == 4);

		// Turn allele frequencies into AD format; order will need to match REF+ALT ordering of nucleotides
		std::vector<int32_t> ad_counts(num_nodes*refalt_count, hts::bcf::int32_missing);
		std::vector<int32_t> ad_info(type_info.width, 0);
        size_t pos = library_start * refalt_count;

		for(size_t u = 0; u < read_depths.size(); ++u) {
            if(type_info.reference == 4) {
                ad_counts[pos++] = 0;
            }
			for(size_t k = 0; k < type_info.width; ++k) {
                int count = read_depths[u][(int)type_info.indexes[k]];
                ad_counts[pos++] = count;
				ad_info[k+(type_info.reference == 4)] += count;
			}
		}

		// Calculate ADR and ADF Tags
		std::vector<int32_t> adf_counts(num_nodes * refalt_count, hts::bcf::int32_missing);
		std::vector<int32_t> adr_counts(num_nodes * refalt_count, hts::bcf::int32_missing);
		std::vector<int32_t> adf_info(refalt_count, 0);
		std::vector<int32_t> adr_info(refalt_count, 0);

		std::vector<int> qual_ref, qual_alt, pos_ref, pos_alt, base_ref, base_alt;
		qual_ref.reserve(depth_stats.dp);
		qual_alt.reserve(depth_stats.dp);
		pos_ref.reserve(depth_stats.dp);
		pos_alt.reserve(depth_stats.dp);
		base_ref.reserve(depth_stats.dp);
		base_alt.reserve(depth_stats.dp);
		double rms_mq = 0.0;

        // zero out the counts for libraries
		for(size_t u = library_start * refalt_count; u < num_nodes * refalt_count; ++u) {
			adf_counts[u] = 0;
			adr_counts[u] = 0;
		}

        int acgt_to_refalt_allele[4] = {-1,-1,-1,-1};
        for(int i=0;i<type_info.width;++i) {
            acgt_to_refalt_allele[(int)type_info.indexes[i]] = i + (type_info.reference == 4);
        }

		for(size_t u = 0; u < data.size(); ++u) {
			const size_t pos = (library_start + u) * refalt_count;
			for(auto && r : data[u]) {
				if(filter_read(r)) {
					continue;
				}
				const size_t base_refalt = acgt_to_refalt_allele[seq::base_index(r.base())];
				assert(base_refalt != -1);
				const size_t base_pos = pos + base_refalt;
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

        record.samples("AD", ad_counts);
		record.samples("ADF", adf_counts);
		record.samples("ADR", adr_counts);

        record.info("AD", ad_info);
		record.info("ADF", adf_info);
		record.info("ADR", adr_info);
		record.info("MQ", static_cast<float>(rms_mq));

        int a11 = adf_info[0];
        int a21 = adr_info[0];
        int a12 = 0, a22 = 0;
        for(int k = 1; k < refalt_count; ++k) {
            a12 += adf_info[k];
            a22 += adr_info[k];
        }
        if(a12+a22 > 0) {
            // Fisher Exact Test for strand bias
            double fs_info = dng::stats::fisher_exact_test(a11, a12, a21, a22);

            double mq_info = dng::stats::ad_two_sample_test(qual_ref, qual_alt);
            double rp_info = dng::stats::ad_two_sample_test(pos_ref, pos_alt);
            double bq_info = dng::stats::ad_two_sample_test(base_ref, base_alt);

            record.info("FS", static_cast<float>(phred(fs_info)));
            record.info("MQTa", static_cast<float>(mq_info));
            record.info("RPTa", static_cast<float>(rp_info));
            record.info("BQTa", static_cast<float>(bq_info));
        }

		record.target(h->target_name[contig]);
		record.position(position);
		vcfout.WriteRecord(record);
		record.Clear();
	});

    return EXIT_SUCCESS;
}

// Process vcf, bcf input data
int process_bcf(task::Call::argument_type &arg) {
    // Read input data
    auto mpileup = io::BcfPileup::open_and_setup(arg);

    auto relationship_graph = create_relationship_graph(arg, &mpileup);

    // Open Output
    auto vcfout = open_vcf_output(arg, mpileup, relationship_graph);

    // Record for each output
    auto record = vcfout.InitVariant();

    // Read header from first file
    const bcf_hdr_t *header = mpileup.reader().header(0); // TODO: fixthis
    const int num_libs = mpileup.num_libraries();

    CallMutations do_call(arg.min_prob, relationship_graph, get_model_parameters(arg));

    // Calculated stats
    CallMutations::stats_t stats;

    // Parameters used by site calculation function
    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;

    // allocate space for ad. bcf_get_format_int32 uses realloc internally
    int n_ad_capacity = num_libs*5;
    auto ad = hts::bcf::make_buffer<int>(n_ad_capacity);

    const double min_prob = arg.min_prob;
    // run calculation based on the depths at each site.
    mpileup([&](const decltype(mpileup)::data_type & rec) {
        // Read all the Allele Depths for every sample into an AD array
        const int n_ad = hts::bcf::get_format_int32(header, rec, "AD", &ad, &n_ad_capacity);
        if(n_ad <= 0) {
            // AD tag is missing, so we do nothing at this time
            // TODO: support using calculated genotype likelihoods
            return;
        }
        assert(n_ad % num_libs == 0);
        const int n_sz = n_ad / num_libs;

        // replace "missing" and "end" values with 0
        boost::replace_if(boost::make_iterator_range(ad.get(), ad.get()+n_ad),
            [](int a) { return a < 0; }, 0 );

        pileup::allele_depths_ref_t read_depths(ad.get(), make_array(num_libs,n_sz));

        if(!do_call(read_depths, n_sz-1, true, &stats)) {
            return;
        }
        const bool has_single_mut = ((stats.mu1p / stats.mup) >= min_prob);
        
        // Set alleles
        record.alleles(static_cast<const char**>((void*)rec->d.allele),
            (int)rec->n_allele);

        // Measure total depth and sort nucleotides in descending order
        pileup::stats_t depth_stats;
        pileup::calculate_stats(read_depths, &depth_stats);

        add_stats_to_output(stats, depth_stats, has_single_mut, relationship_graph, do_call.work(), &record);

        // Turn allele frequencies into AD format; order will need to match REF+ALT ordering of nucleotides
        std::vector<int32_t> ad_counts(num_nodes*n_sz, hts::bcf::int32_missing);
        std::vector<int32_t> ad_info(n_sz, 0);
        size_t pos = library_start * n_sz;

        for(size_t u = 0; u < read_depths.size(); ++u) {
            for(size_t k = 0; k < read_depths[u].size(); ++k) {
                int count = read_depths[u][k];
                ad_counts[pos++] = count;
                ad_info[k] += count;
            }
        }

        record.samples("AD", ad_counts);

        record.info("AD", ad_info);

        // Calculate target position and fetch sequence name
        int contig =  rec->rid;
        int position = rec->pos;

        record.target(mpileup.contigs()[contig].name.c_str());
        record.position(position);
        vcfout.WriteRecord(record);
        record.Clear();
    });
    return EXIT_SUCCESS;
}

int process_ad(task::Call::argument_type &arg) {
#if 0
	// Parse the pedigree file
    Pedigree ped = io::parse_ped(arg.ped);

    // Parse Nucleotide Frequencies
    std::array<double, 4> freqs = utility::parse_nuc_freqs(arg.nuc_freqs);

    // replace arg.region with the contents of a file if needed
    //auto region_ext = io::at_slurp(arg.region);
    if(!arg.region.empty()) {
        throw std::invalid_argument("--region not supported when processing ad/tad file.");
    }

    // Open input files
    if(arg.input.size() != 1) {
        throw std::runtime_error("Can only process one ad/tad file at a time.");
    }
    using AdPileup = dng::io::AdPileup;
    AdPileup mpileup{arg.input[0], std::ios_base::in};

    // Construct peeling algorithm from parameters and pedigree information
    RelationshipGraph relationship_graph;
    if (!relationship_graph.Construct(ped, mpileup.libraries(), inheritance_model(arg.model),
                                      arg.mu, arg.mu_somatic, arg.mu_library,
                                      arg.normalize_somatic_trees)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop relationship_graph.");
    }

    // Select libraries in the input that are used in the pedigree
    mpileup.SelectLibraries(relationship_graph.library_names());

    // Begin writing VCF header
    auto out_file = vcf_get_output_mode(arg);
    hts::bcf::File vcfout(out_file.first.c_str(), out_file.second.c_str());
    vcf_add_header_text(vcfout, arg);

    for(auto && contig : mpileup.contigs()){
        vcfout.AddContig(contig.name.c_str(), contig.length); // Add contigs to header
    }
    for(auto && line : relationship_graph.BCFHeaderLines()) {
        vcfout.AddHeaderMetadata(line.c_str());  // Add pedigree info
    }
    for(auto && str : relationship_graph.labels()) {
        vcfout.AddSample(str.c_str()); // Add genotype columns
    }
    vcfout.WriteHeader();

    // Record for each output
    auto record = vcfout.InitVariant();

    // Construct Calling Object
    const double min_prob = arg.min_prob;
    CallMutations do_call(min_prob, relationship_graph, {arg.theta, freqs,
            arg.ref_weight, {arg.lib_overdisp, arg.lib_error, arg.lib_bias} });

    // Calculated stats
    CallMutations::stats_t stats;
  
    // Parameters used by site calculation function
    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;

    mpileup([&](const AdPileup::data_type & data) {
        if(!do_call(data, &stats)) {
            return;
        }
        const bool has_single_mut = ((stats.mu1p / stats.mup) >= min_prob);

        // Measure total depth and sort nucleotides in descending order
        pileup::stats_t depth_stats;
        pileup::calculate_stats(data, &depth_stats);

        add_stats_to_output(stats, depth_stats, has_single_mut, relationship_graph, do_call.work(), &record);

        const int color = depth_stats.color;
        const auto & type_info_gt = dng::pileup::AlleleDepths::type_info_gt_table[color];
        const auto & type_info = dng::pileup::AlleleDepths::type_info_table[color];
        const int refalt_count = type_info.width + (type_info.reference == 4);

        // Turn allele frequencies into AD format; order will need to match REF+ALT ordering of nucleotides
        std::vector<int32_t> ad_counts(num_nodes*refalt_count, hts::bcf::int32_missing);
        std::vector<int32_t> ad_info(type_info.width, 0);
        size_t pos = library_start * refalt_count;

        for(size_t u = 0; u < data.num_libraries(); ++u) {
            if(type_info.reference == 4) {
                ad_counts[pos++] = 0;
            }
            for(size_t k = 0; k < data.num_nucleotides(); ++k) {
                int count = data(u,k);
                ad_counts[pos++] = count;
                ad_info[k+(type_info.reference == 4)] += count;
            }
        }

        record.samples("AD", ad_counts);

        record.info("AD", ad_info);

        // Calculate target position and fetch sequence name
        int contig = utility::location_to_contig(data.location());
        int position = utility::location_to_position(data.location());

        record.target(mpileup.contigs()[contig].name.c_str());
        record.position(position);
        vcfout.WriteRecord(record);
        record.Clear();

    });
#endif
    return EXIT_SUCCESS;
}

void append_genotype(const hts::bcf::Variant& rec, int gt, int ploidy, std::string* str) {
    assert(str != nullptr);
    assert(ploidy == 1 || ploidy == 2);
    assert(0 <= gt);

    using namespace hts::bcf;

    if(ploidy == 1) {
        str->append(rec.allele(gt));
    } else {
        auto ab = alleles_from_genotype(gt);
        str->append(rec.allele(ab.first));
        str->append("/");
        str->append(rec.allele(ab.second));
    }
}

void add_stats_to_output(const CallMutations::stats_t& call_stats, const pileup::stats_t& depth_stats,
    bool has_single_mut,
    const RelationshipGraph &graph,    
    const peel::workspace_t &work,
    hts::bcf::Variant *record)
{
    using namespace hts::bcf;

    assert(record != nullptr);
    const size_t num_nodes = work.num_nodes;
    const size_t num_libraries = work.library_nodes.second-work.library_nodes.first;

    record->info("MUP", static_cast<float>(call_stats.mup));
    record->info("LLD", static_cast<float>(call_stats.lld));
    record->info("LLS", static_cast<float>(call_stats.lld-depth_stats.log_null));
    record->info("MUX", static_cast<float>(call_stats.mux));
    record->info("MU1P", static_cast<float>(call_stats.mu1p));

    // Output statistics that are only informative if there is a signal of 1 mutation.
    if(has_single_mut) {
        std::string dnt;
        size_t pos = call_stats.dnl;
        int sz = record->num_alleles();

        if(graph.transitions()[pos].type == dng::RelationshipGraph::TransitionType::Trio) {
            assert(work.ploidies[pos] == 2);
            size_t dad = graph.transition(pos).parent1;
            size_t dad_ploidy = work.ploidies[dad];
            size_t mom = graph.transition(pos).parent2;
            size_t mom_ploidy = work.ploidies[mom];

            size_t width = (mom_ploidy == 2) ? sz*(sz+1)/2 : sz;
            dad = call_stats.dnt_row / width;
            mom = call_stats.dnt_row % width;

            append_genotype(*record, dad, dad_ploidy, &dnt);
            dnt += "x";
            append_genotype(*record, mom, mom_ploidy, &dnt);
            dnt += "->";
            append_genotype(*record, call_stats.dnt_col, work.ploidies[pos], &dnt);
        } else {
            size_t par = graph.transition(pos).parent1;
            append_genotype(*record, par, work.ploidies[par], &dnt);
            dnt += "->";
            append_genotype(*record, call_stats.dnt_col, work.ploidies[pos], &dnt);
        }
        record->info("DNT", dnt);

        record->info("DNL", graph.label(pos));
        record->info("DNQ", call_stats.dnq);
    }

    record->info("DP", depth_stats.dp);

    std::vector<float> float_vector;
    std::vector<int32_t> int32_vector;

    int32_vector.assign(2*num_nodes, hts::bcf::int32_missing);

    assert(call_stats.best_genotypes.size() == num_nodes);
    for(size_t i=0;i<num_nodes;++i) {
        assert(work.ploidies[i] == 1 || work.ploidies[i] == 2);
        auto best = call_stats.best_genotypes[i];
        if(work.ploidies[i] == 2) {
            auto ab = alleles_from_genotype(best);
            int32_vector[2*i] = encode_allele_unphased(ab.first);
            int32_vector[2*i+1] = encode_allele_unphased(ab.second);
        } else {
            int32_vector[2*i] = encode_allele_unphased(best);
            int32_vector[2*i+1] = int32_vector_end;
        }
    }
    record->sample_genotypes(int32_vector);
    record->samples("GQ", call_stats.genotype_qualities);

    const size_t num_alleles = record->num_alleles();
    const size_t gt_width = num_alleles*(num_alleles+1)/2;
    const size_t gt_count = gt_width*num_nodes;
    float_vector.assign(gt_count, float_missing);

    for(size_t i=0,k=0;i<num_nodes;++i) {
        if(work.ploidies[i] == 2) {
            for(size_t j=0;j<gt_width;++j) {
                float_vector[k++] = call_stats.posterior_probabilities[i][j];
            }
        } else if(work.ploidies[i] == 1) {
            size_t j;
            for(j=0;j<num_alleles;++j) {
                float_vector[k++] = call_stats.posterior_probabilities[i][j];
            }
            for(;j<gt_width;++j) {
                float_vector[k++] = float_vector_end;
            }
        }
    }
    record->samples("GP", float_vector);

    float_vector.assign(call_stats.node_mup.begin(), call_stats.node_mup.end());
    record->samples("MUP", float_vector);

    if(has_single_mut) {
        float_vector.assign(call_stats.node_mu1p.begin(), call_stats.node_mu1p.end());
        record->samples("MU1P", float_vector);
    }

    float_vector.assign(gt_count, float_missing);
    for(size_t i=0,k=work.library_nodes.first*gt_width;i<num_libraries;++i) {
        if(work.ploidies[i] == 2) {
            for(size_t j=0;j<gt_width;++j) {
                float_vector[k++] = call_stats.genotype_likelihoods[i][j];
            }
        } else if(work.ploidies[i] == 1) {
            size_t j;
            for(j=0;j<num_alleles;++j) {
                float_vector[k++] = call_stats.genotype_likelihoods[i][j];
            }
            for(;j<gt_width;++j) {
                float_vector[k++] = float_vector_end;
            }

        }        
    }    
    record->samples("GL", float_vector);

    int32_vector.assign(num_nodes, int32_missing);
    for(size_t i=0,k=work.library_nodes.first;i<num_libraries;++i) {
        int32_vector[k++] = depth_stats.node_dp[i];
    }    
    record->samples("DP", int32_vector);
}
