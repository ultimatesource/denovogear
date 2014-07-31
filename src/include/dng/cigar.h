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

#include <htslib/sam.h>

namespace dng { namespace cigar {

inline std::size_t query_length(std::size_t sz, const uint32_t *cigar) {
	std::size_t r = 0;
	for(std::size_t k = 0; k < sz; ++k) {
		r += bam_cigar_oplen(cigar[k])*
			(bam_cigar_type(bam_cigar_op(cigar[k]))&1);
	}
	return r;
}

inline std::size_t target_length(std::size_t sz, const uint32_t *cigar) {
	std::size_t r = 0;
	for(std::size_t k = 0; k < sz; ++k) {
		r += bam_cigar_oplen(cigar[k])*
			((bam_cigar_type(bam_cigar_op(cigar[k]))&2)/2);
	}
	return r;
}

inline uint64_t target_to_query(uint64_t target, uint64_t beg,
	std::size_t sz, const uint32_t *cigar)
{
	assert(target >= beg);
	uint64_t off = target-beg;
	uint64_t pos = 0;
	// scan through cigar string
	for(std::size_t k = 0; k < sz; ++k) {
		// how much of the target does the next cigar op consume?
		uint64_t t = bam_cigar_oplen(cigar[k])*
			((bam_cigar_type(bam_cigar_op(cigar[k]))&2)/2);
		// if we will go past the target, adjust the length to match the off
		if(t > off) {
			pos += off*(bam_cigar_type(bam_cigar_op(cigar[k]))&1);
			break;
		}
		// update the query position
		pos += bam_cigar_oplen(cigar[k])*
			(bam_cigar_type(bam_cigar_op(cigar[k]))&1);
		// consume characters in the offset
		off -= t;
	}
	return pos;
}

}} //namespace dng::cigar

#endif // DNG_CIGAR_H
