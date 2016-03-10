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

// 0-based region specification
struct region_t {
     region_t(location_t b, location_t e) : beg{b}, end{e} { }
     region_t(int tid, int b, int e) : beg{utility::make_location(tid,b)}, end{utility::make_location(tid,e)} { }
     region_t(hts::bam::region_t b) : region_t(b.tid, b.beg, b.end) { }
     
     location_t beg;
     location_t end;
};

// Parse dng region string into bam regions
hts::bam::regions_t bam_parse(const std::string &text, const hts::bam::File &file);

namespace detail {
struct raw_parsed_region_t {
    std::string target;
    int beg; // 1-based
    int end; // 1-based
};

typedef std::vector<raw_parsed_region_t> raw_parsed_regions_t;
std::pair<raw_parsed_regions_t,bool> parse_regions(const std::string &text);
} // namespace dng::regions::detail

} // namespace dng::regions
} // namespace dng

#endif // DNG_REGIONS_H
