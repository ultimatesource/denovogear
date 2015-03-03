/*
 * Copyright (c) 2014-2015 Reed A. Cartwright <reed@cartwrig.ht>
 * Copyright (c) 2015 Kael Dai <kdai1@asu.edu>
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

#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm/replace.hpp>

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

#include <htslib/faidx.h>

#include "version.h"

using namespace dng::task;
using namespace dng;

// Helper function that mimics boost::istream_range
template<class Elem, class Traits> inline
    boost::iterator_range<std::istreambuf_iterator<Elem, Traits> >
istreambuf_range(std::basic_istream<Elem, Traits>& in)
{
    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
        std::istreambuf_iterator<Elem, Traits>(in),
        std::istreambuf_iterator<Elem, Traits>());
}

// Helper function to determines if output should be bcf file, vcf file, or stdout. Also
// parses filename "bcf:<file>" --> "<file>"
std::pair<std::string,std::string> vcf_get_output_mode(Call::argument_type &arg) {
	using boost::algorithm::iequals;

	if(arg.output.empty() || arg.output == "-")
		return {"-","w"};
	auto ret = hts::extra::extract_file_type(arg.output);
	if(iequals(ret.first,"bcf")) {
		return {ret.second,"wb"};
	} else if(iequals(ret.first,"vcf")) {
		return {ret.second,"w"};
	} else {
		throw std::runtime_error("Unknown file format '" + ret.second + "' for output '" + arg.output + "'.");
	}
	return {};
}

// Helper function for writing the necessary 
void vcf_add_header_text(hts::bcf::File &vcfout, Call::argument_type &arg) {
#define XM(lname, sname, desc, type, def) \
	vcfout.AddHeaderMetadata(XS(lname), arg.XV(lname));
#	include <dng/task/call.xmh>
#undef XM	

	// Add the available tags for INFO, FILTER, and FORMAT fields
	// TODO: The commented lines are standard VCF fields that may be worth adding to dng output
	vcfout.AddHeaderMetadata("##INFO=<ID=LL,Number=1,Type=Float,Description=\"Log likelihood\">");
	vcfout.AddHeaderMetadata("##INFO=<ID=PMUT,Number=1,Type=Float,Description=\"Probability of mutation\">");
	//vcfout.AddHeaderField("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
	//vcfout.AddHeaderField("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
	//vcfout.AddHeaderField("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
	//vcfout.AddHeaderField("##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
	vcfout.AddHeaderMetadata("##FILTER=<ID=PASS,Description=\"All filters passed\">");
	//vcfout.AddHeaderField("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	//vcfout.AddHeaderField("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
	//vcfout.AddHeaderField("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
	//vcfout.AddHeaderField("##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");

	// AD defined http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk
	vcfout.AddHeaderMetadata("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
}

void vcf_add_record(hts::bcf::File &vcfout, const char *chrom, int pos, const char ref,
	double ll, double pmut, const std::vector<dng::depth5_t> &read_depths)
{
	vcfout.SetID(chrom, pos, nullptr);
	// vcfout.setQuality() - Need phred scalled measure of the site samples
	vcfout.SetFilter("PASS"); // Currently all sites that make it this far pass the threshold

	vcfout.UpdateInfo("LL", static_cast<float>(ll));
	vcfout.UpdateInfo("PMUT", static_cast<float>(pmut));

	// Based on all the samples determine what nucleotides show up and 
	// the order they will appear in the REF and ALT field
	std::vector<uint16_t> allele_order; // List of nucleotides as they appear in the REF and ALT fields
	std::string allele_order_str; // alt_order in string format, used for SetAlleles

	std::size_t ref_index = seq::char_index(ref);
	allele_order.push_back(ref_index);
	allele_order_str = ref;
	for(std::size_t index = 1; index < 4; index++) {
		// iterate through the three remaining NTs, if the NT exists in one of the sample add it to alt_order
		int alt_allele_index = (ref_index + index)%4;
		for(std::size_t sample = 0; sample < read_depths.size(); sample++) {
			if(read_depths[sample].counts[alt_allele_index] != 0) {
				allele_order.push_back(alt_allele_index);
				allele_order_str += std::string(",") + seq::indexed_char(alt_allele_index);
				break;
			}
		}
	}

	// Update REF, ALT fields
	vcfout.SetAlleles(allele_order_str);

	// Turn allele frequencies into AD format; order will need to match REF+ALT ordering of nucleotides
	std::vector<int32_t> gtcounts;
	for(std::size_t sample = 0; sample < read_depths.size(); sample++) {
		for(std::size_t nt = 0; nt < allele_order.size(); nt++) {
			size_t allele_index = allele_order[nt];
			gtcounts.push_back(read_depths[sample].counts[allele_index]);
		}
	}
	vcfout.UpdateSamples("AD", gtcounts);

	vcfout.WriteRecord();
}


// The main loop for dng-call application
// argument_type arg holds the processed command line arguments
int Call::operator()(Call::argument_type &arg) {
	using namespace std;
	
	// Parse pedigree from file	
	dng::io::Pedigree ped;

	if(!arg.ped.empty()) {
		ifstream ped_file(arg.ped);
		if(!ped_file.is_open()) {
			throw std::runtime_error(
				"unable to open pedigree file '" + arg.ped + "'."
			);
		}
		ped.Parse(istreambuf_range(ped_file));
	} else {
		throw std::runtime_error("pedigree file was not specified.");
	}
	
	// Open Reference
	faidx_t * fai = nullptr;
	if(!arg.fasta.empty()) {
		fai = fai_load(arg.fasta.c_str());
		if(fai == nullptr)
			throw std::runtime_error("unable to open faidx-indexed reference file '"
				+ arg.fasta + "'.");
	}
	
	// If arg.files is specified, read input files from list
	if(!arg.sam_files.empty()) {
		ParsedList list(arg.sam_files.c_str(), ParsedList::kFile);
		arg.input.clear(); // TODO: Throw error/warning if both are specified?
		for(std::size_t i=0; i < list.Size(); ++i)
			arg.input.emplace_back(list[i]);
	}

	// Open up the input sam files
	vector<hts::bam::File> indata;
	for(auto str : arg.input) {
		indata.emplace_back(str.c_str(), "r", arg.region.c_str(), arg.fasta.c_str(),
			arg.min_mapqual, arg.min_qlen);
	}

	// Read the header form the first file
	const bam_hdr_t *h = indata[0].header();
	
	// Construct read groups from the input data
	dng::ReadGroups rgs(indata);

	// Parse Nucleotide Frequencies
	std::array<double, 4> freqs;
	// TODO: read directly into freqs????  This will need a wrapper that provides an "insert" function.
	// TODO: include the size into the pattern, but this makes it harder to catch the second error.
	{
		auto f = util::parse_double_list(arg.nuc_freqs,',',4);
		if(!f.second ) {
			throw std::runtime_error("Unable to parse nuc-freq option. "
				"It must be a comma separated list of floating-point numbers.");
		}
		if(f.first.size() != 4) {
			throw std::runtime_error("Wrong number of values passed to nuc-freq. "
				"Expected 4; found " + std::to_string(f.first.size()) + ".");
		}
		std::copy(f.first.begin(),f.first.end(),&freqs[0]);
	}

	// Construct peeling algorithm from parameters and pedigree information
	dng::Pedigree peeler;
	peeler.Initialize({arg.theta, arg.mu, arg.mu_somatic, arg.mu_library, arg.ref_weight, freqs});
	if(!peeler.Construct(ped,rgs)) {
		throw std::runtime_error("Unable to construct peeler for pedigree; "
			"possible non-zero-loop pedigree.");
	}
	
	// Begin Pileup Phase
	dng::MPileup mpileup(rgs.groups());

#ifdef DEBUG_STDOUT
	// Old method of printing out results, leaving in for now mainly for testing purposes.
	// TODO: remove once vcf output has been throughly tested 
	cout << "Contig\tPos\tRef\tLL\tPmut";
	for(std::string str : rgs.libraries()) {
		cout << '\t' << boost::replace(str, '\t', '.');
	}
	cout << endl;
#else

	// Write VCF header
	auto out = vcf_get_output_mode(arg);
	hts::bcf::File vcfout(out.first.c_str(), out.second.c_str(), PACKAGE_STRING);
	vcf_add_header_text(vcfout, arg);
	
	// Add each genotype/sample column then save the header
	for(std::string str : rgs.libraries()) {
		std::string sample_name = boost::replace(str, '\t', '.');
		vcfout.AddSample(sample_name.c_str());
	}
	vcfout.WriteHeader();
#endif
	// information to hold reference
	char *ref = nullptr;
	int ref_sz = 0;
	int ref_target_id = -1;
	
	// quality thresholds 
	int min_qual = arg.min_basequal;
	double min_prob = arg.min_prob;
	
	// Model genotype likelihoods as a mixuture of two dirichlet multinomials
	// TODO: control these with parameters
	genotype::DirichletMultinomialMixture genotype_likelihood(
			{0.9,0.001,0.001,1.05}, {0.1,0.01,0.01,1.1});
	

	// Preform pileup on data by passing a lambda-function to the mpileup () operator
	std::vector<depth5_t> read_depths(rgs.libraries().size(),{0,0});
	const char gts[10][3] = {"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};

	mpileup(indata, [&](const dng::MPileup::data_type &data, uint64_t loc)
	{
		// Calculate target position and fetch sequence name
		int target_id = location_to_target(loc);
		int position = location_to_position(loc);
		if(target_id != ref_target_id && fai != nullptr) {
			if(ref != nullptr)
				free(ref);
			ref = faidx_fetch_seq(fai, h->target_name[target_id],
				0, 0x7fffffff, &ref_sz);
			ref_target_id = target_id;
		}

		// Calculate reference base
		char ref_base = (ref && 0 <= position && position < ref_sz) ?
			ref[position] : 'N';
		int ref_index = seq::char_index(ref_base);

		// reset all depth counters
		read_depths.assign(read_depths.size(),{0,0});
		
		// pileup on read counts
		// TODO: handle overflow?
		for(std::size_t u=0;u<data.size();++u) {
			for(auto &r : data[u]) {
				if(r.is_missing || r.qual.first[r.pos] < arg.min_basequal)
					continue;
				read_depths[rgs.library_from_id(u)].counts[
					seq::base_index(r.aln.seq_at(r.pos))] += 1;
			}
		}
		// calculate genotype likelihoods and store in the lower library vector
		double scale = 0.0, stemp;
		for(std::size_t u=0;u<read_depths.size();++u) {
			std::tie(peeler.library_lower(u),stemp) = genotype_likelihood({read_depths[u].key},ref_index);
			scale += stemp;
		}

		// Calculate probabilities
		double d = peeler.CalculateLogLikelihood(ref_index)+scale;
		double p = peeler.CalculateMutProbability(ref_index);
		
		// Skip this site if it does not meet lower probability threshold
		if(p < min_prob)
			return;

		//vcfo.addRecord(h->target_name[target_id], position+1, ref_base, d, p, read_depths);
#ifdef DEBUG_STDOUT 
		// Print position and read information
		// TODO: turn this into VCF format
		cout << h->target_name[target_id]
		     << '\t' << position+1;
		cout << '\t' << ref_base;
		cout << '\t' << d << '\t' << p;

		for(std::size_t u=0;u<read_depths.size();++u) {
			cout << '\t' << read_depths[u].counts[0]
		         << ','  << read_depths[u].counts[1]
		         << ','  << read_depths[u].counts[2]
		         << ','  << read_depths[u].counts[3]		         
		         //<< "/[" << peeler.library_lower(u).transpose() << ']';
		         ;
		}

		cout << endl;
#else
		vcf_add_record(vcfout, h->target_name[target_id], position+1, ref_base, d, p, read_depths);
#endif
		return;
	});
		
	// Cleanup
	// TODO: RAII support these things
	if (ref != nullptr)
		free(ref);
	if (fai != nullptr)
		fai_destroy(fai);
	
	return EXIT_SUCCESS;
}
