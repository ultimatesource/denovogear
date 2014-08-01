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

#include <boost/range/iterator_range.hpp>

#include <dng/task/call.h>
#include <dng/pedigree.h>
#include <dng/pedigree_peeler.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>

#include <htslib/faidx.h>

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

int Call::operator()(Call::argument_type &arg) {
	using namespace std;
	
	dng::Pedigree pedigree;

	// Parse pedigree from file	
	if(!arg.ped.empty()) {
		ifstream ped_file(arg.ped);
		if(!ped_file.is_open()) {
			throw std::runtime_error(
				"unable to open pedigree file '" + arg.ped + "'."
			);
		}
		pedigree.Parse(istreambuf_range(ped_file));
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
	vector<fileio::SamFile> indata;
		
	for(auto str : arg.input) {
		indata.emplace_back(str.c_str(), arg.region.c_str(), arg.fasta.c_str(),
			arg.min_mapqual, arg.min_qlen);
	}
	const bam_hdr_t *h = indata[0].header();
	
	dng::ReadGroups rgs(indata);

	dng::MPileup mpileup(rgs.groups());

	cerr << "Contig\tPos\tRef";
	for(auto a : rgs.groups()) {
		cerr << '\t' << a;
	}
	cerr << endl;
	
	// information to hold reference
	char *ref = nullptr;
	int ref_sz = 0;
	int ref_target_id = -1;
	
	mpileup(indata, [h,fai,&ref,&ref_sz,&ref_target_id](const dng::MPileup::data_type &data, uint64_t loc)
	{
		int target_id = location_to_target(loc);
		int position = location_to_position(loc);
		if(target_id != ref_target_id && fai != nullptr) {
			if(ref != nullptr)
				free(ref);
			ref = faidx_fetch_seq(fai, h->target_name[target_id],
				0, 0x7fffffff, &ref_sz);
			ref_target_id = target_id;
		}
		char ref_base = (ref && 0 <= position && position < ref_sz) ?
			ref[position] : 'N';
		
		cerr << h->target_name[target_id]
		     << "\t" << position+1
		     << "\t" << ref_base;
		for(auto &d : data) {
			cerr << "\t" << d.size();
		}
		cout << "\n";
		return;
	});
	
	/*
	dng::PedigreePeeler peeler;
	peeler.Initialize(arg.theta, arg.mu);
	peeler.Construct(pedigree);
	
	dng::PedigreePeeler::IndividualBuffer buf;
	buf.resize(7, dng::PedigreePeeler::Vector10d::Zero());
	buf[1][1] = 1.0;
	buf[2][1] = 1.0;
	buf[3][1] = 1.0;
	buf[4][0] = 1.0;
	buf[5][0] = 1.0;
	buf[6][0] = 1.0;
	
	double d = peeler.CalculateLogLikelihood(buf);
	cout << exp(d) << endl;
	*/
	if (ref != nullptr)
		free(ref);
	if (fai != nullptr)
		fai_destroy(fai);
	
	
	return EXIT_SUCCESS;
}
