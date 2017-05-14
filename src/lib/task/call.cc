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

#include <boost/algorithm/string.hpp>

#include <dng/task/call.h>
#include <dng/relationship_graph.h>
#include <dng/fileio.h>
#include <dng/read_group.h>
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
#include <dng/io/ped.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>

#include "version.h"

using namespace dng;

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

// Helper function for writing the vcf header information
void vcf_add_header_text(hts::bcf::File &vcfout, task::Call::argument_type &arg) {
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
	// Check gamma parameter
    if(arg.gamma.size() < 2) {
        throw std::runtime_error("Unable to construct genotype-likelihood model; "
                                 "Gamma needs to be specified at least twice to change model from default.");
    }

    // replace arg.region with the contents of a file if needed
    io::at_slurp(arg.region);

    // Determine the type of input files
    auto it = arg.input.begin();
    FileCat mode = utility::input_category(*it, FileCat::Sequence|FileCat::Pileup|FileCat::Variant, FileCat::Unknown);
    for(++it; it != arg.input.end(); ++it) {
    	// Make sure different types of input files aren't mixed together
        if(utility::input_category(*it, FileCat::Sequence|FileCat::Pileup|FileCat::Variant, FileCat::Sequence) != mode) {
            throw std::runtime_error("Mixing pileup, sequencing, and variant file types is not supported.");
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
    	throw std::runtime_error("Unknown input data file type.");
    }
    return EXIT_FAILURE;
}

// Processes bam, sam, and cram files.
int process_bam(task::Call::argument_type &arg) {
    // Parse pedigree from file

    Pedigree ped = io::parse_ped(arg.ped);

    // Open Reference
    io::Fasta reference{arg.fasta.c_str()};

    // Parse Nucleotide Frequencies
    std::array<double, 4> freqs = utility::parse_nuc_freqs(arg.nuc_freqs);

    // Fetch the selected regions
    auto region_ext = io::at_slurp(arg.region); // replace arg.region with the contents of a file if needed

    // Open input files
    using BamPileup = dng::io::BamPileup;
    BamPileup mpileup{arg.min_qlen, arg.rgtag};
    for(auto && str : arg.input) {
        // TODO: We put all this logic here to simplify the BamPileup construction
        // TODO: However it might make sense to incorporate it into BamPileup::AddFile
        hts::bam::File input{str.c_str(), "r", arg.fasta.c_str(), arg.min_mapqual, arg.header.c_str()};
        if(!input.is_open()) {
            throw std::runtime_error("Unable to open input file '" + str + "'.");
        }
        // add regions
        // TODO: make this work with region_ext
        if(!arg.region.empty()) {
            if(arg.region.find(".bed") != std::string::npos) {
                input.regions(regions::bam_parse_bed(arg.region, input));
            } else {
                input.regions(regions::bam_parse_region(arg.region, input));
            }
        }
        mpileup.AddFile(std::move(input));
    }

    // Construct peeling algorithm from parameters and pedigree information
    RelationshipGraph relationship_graph;
    if (!relationship_graph.Construct(ped, mpileup.libraries(), inheritance_model(arg.model),
                                      arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop relationship_graph.");
    }

    // Select libraries in the input that are used in the pedigree
    mpileup.SelectLibraries(relationship_graph.library_names());

    // Write VCF output header
    auto out_file = vcf_get_output_mode(arg);
    hts::bcf::File vcfout(out_file.first.c_str(), out_file.second.c_str());
    vcf_add_header_text(vcfout, arg);

    // User the header from the first file to determine the contigs
    for(auto && contig : mpileup.contigs()) {
    	vcfout.AddContig(contig.name.c_str(), contig.length); // Add contigs to header
    }

    for(auto && line : relationship_graph.BCFHeaderLines()) {
    	vcfout.AddHeaderMetadata(line.c_str()); // Add pedigree information
    }

    // Finish Header
    for(auto && str : relationship_graph.labels()) {
        vcfout.AddSample(str.c_str()); // Create samples
    }
    vcfout.WriteHeader();

    // Record for each output
    auto record = vcfout.InitVariant();

    // Construct Calling Object
    const double min_prob = arg.min_prob;
    CallMutations do_call(min_prob, relationship_graph, {arg.theta, freqs,
            arg.ref_weight, arg.gamma[0], arg.gamma[1]});

    // Pileup data
    dng::pileup::RawDepths read_depths(mpileup.num_libraries());

    // Calculated stats
    CallMutations::stats_t stats;
  
    // Parameters used by site calculation function
    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;
    const int min_basequal = arg.min_basequal;
    auto filter_read = [min_basequal](
    BamPileup::data_type::value_type::const_reference r) -> bool {
        return (r.is_missing
        || r.base_qual() < min_basequal
        || seq::base_index(r.base()) >= 4);
    };

    auto h = mpileup.header(); // TODO: Fixme.

    mpileup([&](const BamPileup::data_type & data, utility::location_t loc) {
        // Calculate target position and fetch sequence name
        int contig = utility::location_to_contig(loc);
        int position = utility::location_to_position(loc);

        // Calculate reference base
        assert(0 <= contig && contig < h->n_targets);
        char ref_base = reference.FetchBase(h->target_name[contig],position);
        size_t ref_index = seq::char_index(ref_base);

		// reset all depth counters
        boost::fill(read_depths, pileup::depth_t{});

		// pileup on read counts
		for(std::size_t u = 0; u < data.size(); ++u) {
			for(auto && r : data[u]) {
				if(filter_read(r)) {
					continue;
				}
                std::size_t base = seq::base_index(r.base());
                assert(read_depths[u].counts[ base ] < 65535);

                read_depths[u].counts[ base ] += 1;
			}
		}
		if(!do_call(read_depths, ref_index, &stats)) {
			return;
		}
        const bool has_single_mut = ((stats.mu1p / stats.mup) >= min_prob);

		// Measure total depth and sort nucleotides in descending order
        pileup::stats_t depth_stats;
        pileup::calculate_stats(read_depths, ref_index, &depth_stats);

        add_stats_to_output(stats, depth_stats, has_single_mut, relationship_graph, do_call.work(), &record);

        // Information about the site, based on depths
        const int color = depth_stats.color;
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
                int count = read_depths[u].counts[(int)type_info.indexes[k]];
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

        record.samples("AD", ad_counts);
		record.samples("ADF", adf_counts);
		record.samples("ADR", adr_counts);

        record.info("AD", ad_info);
		record.info("ADF", adf_info);
		record.info("ADR", adr_info);
		record.info("MQ", static_cast<float>(rms_mq));
		record.info("FS", static_cast<float>(phred(fs_info)));
		record.info("MQTa", static_cast<float>(mq_info));
		record.info("RPTa", static_cast<float>(rp_info));
		record.info("BQTa", static_cast<float>(bq_info));

		record.target(h->target_name[contig]);
		record.position(position);
		vcfout.WriteRecord(record);
		record.Clear();
	});

    return EXIT_SUCCESS;
}

int process_ad(task::Call::argument_type &arg) {
	// Parse the pedigree file
    Pedigree ped = io::parse_ped(arg.ped);

    // Parse Nucleotide Frequencies
    std::array<double, 4> freqs = utility::parse_nuc_freqs(arg.nuc_freqs);

    //io::at_slurp(arg.region); // replace arg.region with the contents of a file if needed

    // Open input files
    if(arg.input.size() != 1) {
        throw std::runtime_error("Can only process one ad/tad file at a time.");
    }
    using AdPileup = dng::io::AdPileup;
    AdPileup mpileup{arg.input[0], std::ios_base::in};

    // Construct peeling algorithm from parameters and pedigree information
    RelationshipGraph relationship_graph;
    if (!relationship_graph.Construct(ped, mpileup.libraries(), inheritance_model(arg.model),
                                      arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop relationship_graph.");
    }

    // Select libraries in the input that are used in the pedigree
    std::cerr << "before " << mpileup.libraries().names.size() << "\n";
    std::cerr << relationship_graph.library_names().size() << "\n";
    mpileup.SelectLibraries(relationship_graph.library_names());
    std::cerr << "after " << mpileup.libraries().names.size() << "\n";

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
            arg.ref_weight, arg.gamma[0], arg.gamma[1]});

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

    return EXIT_SUCCESS;
}

// Process vcf, bcf input data
int process_bcf(task::Call::argument_type &arg) {
	// Parse the pedigree file
    Pedigree ped = io::parse_ped(arg.ped);

    // Parse Nucleotide Frequencies
    auto freqs = utility::parse_nuc_freqs(arg.nuc_freqs);

    // Read input data
    if(arg.input.size() > 1) {
    	throw std::runtime_error("dng call can only handle one variant file at a time.");
    }

    using dng::io::BcfPileup;
    BcfPileup mpileup;

    // Fetch the selected regions
    auto region_ext = io::at_slurp(arg.region); // replace arg.region with the contents of a file if needed
    if(!arg.region.empty()) {
        auto ranges = regions::parse_ranges(arg.region);
        if(ranges.second == false) {
            throw std::runtime_error("unable to parse the format of '--region' argument.");
        }
        mpileup.SetRegions(ranges.first);
    }

    if(mpileup.AddFile(arg.input[0].c_str()) == 0) {
        int errnum = mpileup.reader().handle()->errnum;
        throw std::runtime_error(bcf_sr_strerror(errnum));
    }

    // Construct peeling algorithm from parameters and pedigree information
    InheritanceModel model = inheritance_model(arg.model);

    dng::RelationshipGraph relationship_graph;
    if (!relationship_graph.Construct(ped, mpileup.libraries(), model,
                                      arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                 "possible non-zero-loop relationship_graph.");
    }
    mpileup.SelectLibraries(relationship_graph.library_names());

    // Begin writing VCF header
    auto out_file = vcf_get_output_mode(arg);
    hts::bcf::File vcfout(out_file.first.c_str(), out_file.second.c_str());
    vcf_add_header_text(vcfout, arg);

    // Read header from first file
    const bcf_hdr_t *header = mpileup.reader().header(0); // TODO: fixthis
    const int num_libs = mpileup.num_libraries();

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

    // Record for each output
    auto record = vcfout.InitVariant();

    double min_prob = arg.min_prob;
    CallMutations do_call(min_prob, relationship_graph, {arg.theta, freqs,
            arg.ref_weight, arg.gamma[0], arg.gamma[1]});

    // Calculated stats
    CallMutations::stats_t stats;

    // Parameters used by site calculation function
    const size_t num_nodes = relationship_graph.num_nodes();
    const size_t library_start = relationship_graph.library_nodes().first;

    using pileup::AlleleDepths;

    AlleleDepths data;

    // allocate space for ad. bcf_get_format_int32 uses realloc internally
    int n_ad_capacity = num_libs*5;
    std::unique_ptr<int[],decltype(std::free)*> ad{
        reinterpret_cast<int*>(std::malloc(sizeof(int)*n_ad_capacity)), std::free};
    if(!ad) {
        throw std::bad_alloc{};
    }

    // run calculation based on the depths at each site.
    mpileup([&](const BcfPileup::data_type & rec) {
        data.location(rec->rid, rec->pos);

        // Identify site
        std::string allele_str;
        std::vector<char> allele_indexes;
        std::vector<int>  allele_list;
        const int num_alleles = rec->n_allele;
        int first_is_n = 0;
        for(int a = 0; a < num_alleles; ++a) {
            int n = AlleleDepths::MatchAlleles(rec->d.allele[a]);
            if(0 <= n && n <= 3) {
                allele_indexes.push_back(n);
                allele_list.push_back(a);
            } else if(a == 0 && n == 4) {
                first_is_n = 64;
            }
        }
        int color = AlleleDepths::MatchIndexes(allele_indexes);
        assert(color != -1);
        data.resize(color+first_is_n, rec->n_sample);
        // Read all the Allele Depths for every sample into an AD array

        int *pad = ad.get();
        int n_ad = bcf_get_format_int32(header, rec, "AD", &pad, &n_ad_capacity);
        if(n_ad == -4) {
            throw std::bad_alloc{};
        } else if(pad != ad.get()) {
            // update pointer
            ad.release();
            ad.reset(pad);
        }
        if(n_ad <= 0) {
            // AD tag is missing, so we do nothing at this time
            // TODO: support using calculated genotype likelihoods
            return;
        }

        assert(n_ad >= data.data_size());
        for(int i=0;i<num_libs;++i) {
            int offset = i*num_alleles;
            for(int a=0; a<allele_indexes.size();++a) {
                data(i,a) = ad[offset+allele_list[a]];
            }
        }

        if(!do_call(data, &stats)) {
            return;
        }
        const bool has_single_mut = ((stats.mu1p / stats.mup) >= min_prob);

        // Measure total depth and sort nucleotides in descending order
        pileup::stats_t depth_stats;
        pileup::calculate_stats(data, &depth_stats);

        add_stats_to_output(stats, depth_stats, has_single_mut, relationship_graph, do_call.work(), &record);

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

    return EXIT_SUCCESS;
}

std::string genotype_string(int gt, int ploidy, int color);

inline
int reencode_genotype(int gt, int ploidy, int from_color, int to_color) {
    assert(0 <= from_color && from_color < pileup::AlleleDepths::type_info_table_length);
    assert(0 <= to_color && to_color < pileup::AlleleDepths::type_info_table_length);
    const auto & type_info_gt_table = dng::pileup::AlleleDepths::type_info_gt_table;
    const auto & type_info_table = dng::pileup::AlleleDepths::type_info_table;
    return (from_color == to_color) ? gt :
           (ploidy == 2) ? utility::find_position(type_info_gt_table[to_color].indexes,
                               type_info_gt_table[from_color].indexes[gt])
                         : utility::find_position(type_info_table[to_color].indexes,
                               type_info_table[from_color].indexes[gt]);
}

void add_stats_to_output(const CallMutations::stats_t& call_stats, const pileup::stats_t& depth_stats,
    bool has_single_mut,
    const RelationshipGraph &graph,    
    const peel::workspace_t &work,
    hts::bcf::Variant *record)
{
    assert(record != nullptr);
    const int old_color = call_stats.color;
    const int new_color = depth_stats.color;
    const size_t num_nodes = work.num_nodes;
    const size_t num_libraries = work.library_nodes.second-work.library_nodes.first;

    const auto & type_info_gt_table = dng::pileup::AlleleDepths::type_info_gt_table;
    const auto & type_info_table = dng::pileup::AlleleDepths::type_info_table;

    record->alleles(type_info_table[new_color].label_htslib);

    record->info("MUP", static_cast<float>(call_stats.mup));
    record->info("LLD", static_cast<float>(call_stats.lld));
    record->info("LLS", static_cast<float>(call_stats.lld-depth_stats.log_null));
    record->info("MUX", static_cast<float>(call_stats.mux));
    record->info("MU1P", static_cast<float>(call_stats.mu1p));

    // Output statistics that are only informative if there is a signal of 1 mutation.
    if(has_single_mut) {
        std::string dnt;
        size_t pos = call_stats.dnl;
        if(graph.transitions()[pos].type == dng::RelationshipGraph::TransitionType::Trio) {
            assert(work.ploidies[pos] == 2);
            size_t child = reencode_genotype(call_stats.dnt_col, 2, old_color, new_color);

            size_t dad = graph.transition(pos).parent1;
            size_t dad_ploidy = work.ploidies[dad];
            size_t mom = graph.transition(pos).parent2;
            size_t mom_ploidy = work.ploidies[mom];

            size_t width = (mom_ploidy == 2) ? type_info_gt_table[old_color].width : type_info_table[new_color].width;
            dad = call_stats.dnt_row / width;
            mom = call_stats.dnt_row % width;
            dad = reencode_genotype(dad, dad_ploidy, old_color, new_color);
            mom = reencode_genotype(mom, mom_ploidy, old_color, new_color);

            dnt = genotype_string(dad, dad_ploidy, new_color);
            dnt += 'x';
            dnt += genotype_string(mom, mom_ploidy, new_color);
            dnt += '>';
            dnt += genotype_string(child, work.ploidies[pos], new_color);            
        } else {
            size_t par_pos = graph.transition(pos).parent1;
            size_t par   = reencode_genotype(call_stats.dnt_row, work.ploidies[par_pos], old_color, new_color);
            size_t child = reencode_genotype(call_stats.dnt_col, work.ploidies[pos], old_color, new_color);

            dnt = genotype_string(par, work.ploidies[par_pos], new_color);
            dnt += '>';
            dnt += genotype_string(child, work.ploidies[pos], new_color);
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
        auto best = call_stats.best_genotypes[i];
        if(work.ploidies[i] == 2) {
            size_t pos = reencode_genotype(best, 2, old_color, new_color);
            assert(pos < type_info_gt_table[new_color].width);
            int32_vector[2*i] = dng::pileup::AlleleDepths::encoded_alleles_diploid_unphased[pos][0];
            int32_vector[2*i+1] = dng::pileup::AlleleDepths::encoded_alleles_diploid_unphased[pos][1];
        } else if(work.ploidies[i] == 1) {
            size_t pos = reencode_genotype(best, 1, old_color, new_color);
            assert(pos < type_info_table[new_color].width);
            int32_vector[2*i] = dng::pileup::AlleleDepths::encoded_alleles_haploid[pos][0];
            int32_vector[2*i+1] = dng::pileup::AlleleDepths::encoded_alleles_haploid[pos][1];
        } else {
            assert(0); // should not be reached, only ploidies of 1 and 2 are supported
        }
    }
    record->sample_genotypes(int32_vector);
    record->samples("GQ", call_stats.genotype_qualities);

    const size_t gt_width = type_info_gt_table[new_color].width;
    const size_t gt_count = gt_width*num_nodes;
    float_vector.assign(gt_count, hts::bcf::float_missing);
    for(size_t i=0,k=0;i<num_nodes;++i) {
        if(work.ploidies[i] == 2) {
            for(size_t j=0;j<gt_width;++j) {
                size_t pos = reencode_genotype(j, 2, new_color, old_color);
                assert(pos < type_info_gt_table[old_color].width);
                float_vector[k++] = call_stats.posterior_probabilities[i][pos];
            }
        } else if(work.ploidies[i] == 1) {
            size_t j;
            for(j=0;j<type_info_table[new_color].width;++j) {
                size_t pos = reencode_genotype(j, 1, new_color, old_color);
                assert(pos < type_info_table[old_color].width);
                float_vector[k++] = call_stats.posterior_probabilities[i][pos];
            }
            for(;j<gt_width;++j) {
                float_vector[k++] = hts::bcf::float_vector_end;
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

    float_vector.assign(gt_count, hts::bcf::float_missing);
    for(size_t i=0,k=work.library_nodes.first*gt_width;i<num_libraries;++i) {
        if(work.ploidies[i] == 2) {
            for(size_t j=0;j<gt_width;++j) {
                size_t pos = reencode_genotype(j, 2, new_color, old_color);
                assert(pos < type_info_gt_table[old_color].width);
                float_vector[k++] = call_stats.genotype_likelihoods[i][pos];
            }
        } else if(work.ploidies[i] == 1) {
            size_t j;
            for(j=0;j<type_info_table[new_color].width;++j) {
                size_t pos = reencode_genotype(j, 1, new_color, old_color);
                assert(pos < type_info_table[old_color].width);
                float_vector[k++] = call_stats.genotype_likelihoods[i][pos];
            }
            for(;j<gt_width;++j) {
                float_vector[k++] = hts::bcf::float_vector_end;
            }

        }        
    }    
    record->samples("GL", float_vector);

    int32_vector.assign(num_nodes, hts::bcf::int32_missing);
    for(size_t i=0,k=work.library_nodes.first;i<num_libraries;++i) {
        int32_vector[k++] = depth_stats.node_dp[i];
    }    
    record->samples("DP", int32_vector);
}

std::string genotype_string(int gt, int ploidy, int color) {
    using AlleleDepths = dng::pileup::AlleleDepths;
    assert(0 <= color && color < AlleleDepths::type_info_table_length);
    auto & type_info = AlleleDepths::type_info_table[color];
    auto & type_info_gt = AlleleDepths::type_info_gt_table[color];

    std::string out;
    if(ploidy == 2) {
        int a = AlleleDepths::alleles_diploid[gt][0] + (type_info.reference == 4);
        int b = AlleleDepths::alleles_diploid[gt][1] + (type_info.reference == 4);
        assert(0 <= a && a < strlen(type_info.label_upper));
        assert(0 <= b && b < strlen(type_info.label_upper));
        out += type_info.label_upper[a];
        out += type_info.label_upper[b];
    } else if(ploidy == 1) {
        int a = gt + (type_info.reference == 4);
        assert(0 <= a && a < strlen(type_info.label_upper));        
        out += type_info.label_upper[a];
    } else {
        assert(0); // should not reach here
    }
    return out;
}
