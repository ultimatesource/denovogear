/*
 * Copyright (c) 2016 Reed A. Cartwright
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
#ifndef DNG_REGIONS_H
#define DNG_REGIONS_H

#include <utility>
#include <string>

#include <dng/utility.h>
#include <dng/hts/bam.h>

#include <dng/detail/unit_test.h>

namespace dng {
namespace regions {

// Parse a white-spaced separated list of regions
// Formats of a region range is contig:from-to
//   contig: name of contig/fragment/sequence
//   from: beginning of range, 1-based inclusive
//   to: end of range, 1-based inclusive
// Format list:
//   ctg ctg:from-to ctg:position ctg:from- ctg:-to
//   ctg:from+length

struct parsed_range_t {
    std::string target;
    int beg; // 1-based
    int end; // 1-based
};

typedef std::vector<parsed_range_t> parsed_ranges_t;

std::pair<parsed_ranges_t,bool> parse_ranges(const std::string &text);

// 0-based region specification
struct location_range_t {
     location_range_t(location_t b, location_t e) : beg{b}, end{e} { }
     location_range_t(int tid, int b, int e) : beg{utility::make_location(tid,b)}, end{utility::make_location(tid,e)} { }
     location_range_t(hts::bam::region_t b) : location_range_t(b.tid, b.beg, b.end) { }
     
     location_t beg;
     location_t end;
};

// Parse dng region string into bam regions
hts::bam::regions_t bam_parse_region(const std::string &text, const hts::bam::File &file);

// Parse dng region from bed into bam regions
hts::bam::regions_t bam_parse_bed(const std::string &text, const hts::bam::File &file);

} // namespace dng::regions
} // namespace dng

#endif // DNG_REGIONS_H
