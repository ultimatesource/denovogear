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

#pragma once
#ifndef DNG_PILEUP_H
#define DNG_PILEUP_H

#include <vector>
#include <limits>
#include <htslib/sam.h>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/metafunctions.hpp> 

#include <dng/fileio.h>

namespace dng {

typedef void (mpileup_func_t)(int,int,const std::vector<int>&,
	const std::vector<const bam_pileup1_t *>&);

template<typename Input, typename Func>
void mpileup(Input &range, bam_plp_auto_f auto_func, Func func) {
	using namespace std;
	using namespace fileio;
		
	// type erase our callback function
	function<mpileup_func_t> call_back(func);
	
	// bam_mplp_init requires an array of pointers
	std::size_t num_files = boost::size(range);
	vector<typename boost::range_pointer<Input>::type>
		pointers(num_files, nullptr);
	{ std::size_t i=0;
	for(auto a=boost::begin(range); a != boost::end(range); ++a) {
		pointers[i] = &(*a) + i;
		i++;
	}}
	
	// If there is a parsed region, use it.
	// TODO: check to see if the regions are all the same???
	int beg = 0, end = std::numeric_limits<int>::max();
	if(boost::begin(range)->iter() != nullptr) {
		beg = boost::begin(range)->iter()->beg;
		end = boost::begin(range)->iter()->end;
	}
	
	// initialize the pileup iterator
	bam_mplp_t mpileup_iter = bam_mplp_init(num_files, auto_func,
		reinterpret_cast<void**>(pointers.data()));
	
	vector<int> counts(num_files, 0);
	vector<const bam_pileup1_t *> records(num_files, nullptr);
	int target_id, pos;
	
	while(bam_mplp_auto(mpileup_iter, &target_id, &pos, counts.data(),records.data()) > 0) {
		if (pos < beg || pos >= end) continue; // out of range; skip
		//if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
		call_back(target_id, pos, counts, records);
	}
	bam_mplp_destroy(mpileup_iter);
}

} //namespace dng

#endif //DNG_PILEUP_H

