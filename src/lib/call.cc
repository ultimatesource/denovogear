/*
 * Copyright (c) 2014 Reed A. Cartwright
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

#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm/replace.hpp>

#include <boost/spirit/home/x3.hpp>

#include <boost/algorithm/string.hpp>

#include <dng/task/call.h>
#include <dng/pedigree.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/likelihood.h>
#include <dng/seq.h>
#include <dng/hts/vcf.h>

#include <htslib/faidx.h>

using namespace dng::task;
using namespace dng;
namespace x3 = boost::spirit::x3;

// Helper function that mimics boost::istream_range
template<class Elem, class Traits> inline
    boost::iterator_range<std::istreambuf_iterator<Elem, Traits> >
istreambuf_range(std::basic_istream<Elem, Traits>& in)
{
    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
        std::istreambuf_iterator<Elem, Traits>(in),
        std::istreambuf_iterator<Elem, Traits>());
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
	// TODO: turn this into utility function
	{
		std::vector<double> f;
		f.reserve(4);
		x3::ascii::space_type space;
		auto b = arg.nuc_freqs.begin();
		auto e = arg.nuc_freqs.end();
		bool r = x3::phrase_parse(b, e, x3::double_ % ',', space,f);
		if(!r || b != e ) {
			throw std::runtime_error("Unable to parse nuc-freq option. It must be a comma separated list of floating-point numbers.");
		}
		if(f.size() != 4) {
			throw std::runtime_error("Wrong number of values passed to nuc-freq. Expected 4; found " + std::to_string(f.size()) + ".");
		}
		std::copy(f.begin(),f.end(),&freqs[0]);
	}

	// Construct peeling algorithm from parameters and pedigree information
	dng::Pedigree peeler;
	peeler.Initialize({arg.theta, arg.mu, arg.mu_somatic, arg.mu_library, arg.ref_weight, freqs});
	if(!peeler.Construct(ped,rgs)) {
		throw std::runtime_error("Unable to construct peeler for pedigree; possible non-zero-loop pedigree.");
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
	string mode("w");
	if(boost::algorithm::ends_with(boost::to_upper_copy(arg.output), ".BCF") == true)
	  {
	    // if file ends with .bcf write to a bcf file instead of VCF
	    mode += "b";
	  }
	// TODO: Add code to allow for compressed output

	hts::vcf::Vcfostream vcfo(arg.output.c_str(), mode.c_str(), arg.ped.c_str());
	//hts::vcf::Vcfostream vcfo("-", "w", "ceu.ped");
	//vcfo.addHdrMetadata("fasta", arg.fasta.c_str());
	vcfo.addHdrMetadata("region", arg.region.c_str());
	vcfo.addHdrMetadata("min_qlen", arg.min_qlen);
	vcfo.addHdrMetadata("min_mapqual", arg.min_mapqual);
	vcfo.addHdrMetadata("mu", arg.mu);	
	vcfo.addHdrMetadata("mu_somatic", arg.mu_somatic);
	vcfo.addHdrMetadata("mu_library", arg.mu_library);
	vcfo.addHdrMetadata("theta", arg.theta);
	vcfo.addHdrMetadata("ref_weight", arg.ref_weight);
	vcfo.addHdrMetadata("min_basequal", arg.min_basequal);
	
	// Add each genotype column
	for(std::string str : rgs.libraries())
	  {
	    std::string sample_name = boost::replace(str, '\t', '.');
	    vcfo.addSample(sample_name.c_str());
	  }
	vcfo.writeHdr();
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
		vcfo.addRecord(h->target_name[target_id], position+1, ref_base, d, p, read_depths);
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
