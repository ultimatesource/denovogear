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
#ifndef DNG_CIGAR_H
#define DNG_CIGAR_H

#include <cstdlib>
#include <cassert>
#include <utility>

#include <htslib/sam.h>

namespace dng { namespace cigar {

typedef std::pair<std::size_t, const uint32_t *> cigar_t;

inline std::size_t query_length(cigar_t cigar) {
	std::size_t r = 0;
	for(std::size_t k = 0; k < cigar.first; ++k) {
		r += bam_cigar_oplen(cigar.second[k])*
			(bam_cigar_type(bam_cigar_op(cigar.second[k]))&1);
	}
	return r;
}

inline std::size_t target_length(cigar_t cigar) {
	std::size_t r = 0;
	for(std::size_t k = 0; k < cigar.first; ++k) {
		r += bam_cigar_oplen(cigar.second[k])*
			((bam_cigar_type(bam_cigar_op(cigar.second[k]))&2)/2);
	}
	return r;
}


// Position in query is result/2.
// Result&0x1 == 1 if query contants a gap at target position.

inline uint64_t query_pos(uint64_t q) { return q/2; }
inline uint64_t query_del(uint64_t q) { return q&1; }

inline uint64_t target_to_query(uint64_t target, uint64_t beg,
	cigar_t cigar)
{
	if(target < beg)
		return -1;
	uint64_t off = (target-beg)*2;
	uint64_t pos = -1;
	// scan through cigar string
	for(std::size_t k = 0; k < cigar.first; ++k) {
		// how much of the target does the next cigar op consume?
		uint64_t t = bam_cigar_oplen(cigar.second[k])*
			((bam_cigar_type(bam_cigar_op(cigar.second[k]))&2));
		// if we will go past the target, adjust the length to match the off
		if(t > off)
			return pos+(bam_cigar_type(bam_cigar_op(cigar.second[k]))&1)*(off+1);
		// update the query position
		pos += (bam_cigar_oplen(cigar.second[k]))*
			(bam_cigar_type(bam_cigar_op(cigar.second[k]))&1)*2;
		// consume characters in the offset
		off -= t;
	}
	return pos;
}

} // namespace cigar
using cigar::cigar_t;
} //namespace dng

#endif // DNG_CIGAR_H
