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
#include <vector>
#include <queue>
#include <deque>

#include <dng/utility.h>
#include <dng/io/utility.h>

#include <boost/container/flat_map.hpp>
#include <boost/optional.hpp>

#include <dng/detail/unit_test.h>

namespace dng {
namespace regions {

struct contig_t {
    contig_t() = default;
    contig_t(std::string s, int n) : name{std::move(s)}, length{n} { }
    std::string name;
    int length;
};

class ContigIndex {
public:
    int NameToId(const std::string &name) const {
        auto it = map_.find(name);
        return (it != map_.end()) ? it->second : -1; 
    }

    int NameToId(const std::string &name, int hint) const {
        assert(0 <= hint && hint < contigs_.size());
        if(contigs_[hint].name == name) {
            return hint;
        } else if(hint+1 < contigs_.size() && contigs_[hint+1].name == name) {
            return hint+1;
        }
        return NameToId(name);
    }

    const contig_t* IdToContig(int id) const {
        return (0 <= id && id < contigs_.size()) ? &contigs_[id] : nullptr;
    }

    const std::vector<contig_t>& contigs() const { return contigs_; }

    const contig_t& contig(int id) const { return contigs_[id]; }

    template<typename... Args>
    int AddContig(Args&&... args) {
        int pos = contigs_.size();
        // Build the contig early because we need to access its name
        contigs_.emplace_back(std::forward<Args>(args)...);
        // Try to insert it into the name map
        auto ok = map_.emplace(contigs_.back().name, pos);
        if(ok.second) {
            return pos;
        }
        // Reject the contig because the name is a duplicate
        contigs_.pop_back();
        return -1;
    }

private:
    // Keeps contigs sorted in file order.
    std::vector<contig_t> contigs_;

    // Allows fast conversion between contig_name and id
    boost::container::flat_map<std::string, int> map_;
};

// 1-based or 0-based region specification with named contig
struct contig_fragment_t {
    std::string contig_name;
    int beg;
    int end;
};

typedef std::vector<contig_fragment_t> contig_fragments_t;

// 0-based region specification
struct range_t {
     range_t(location_t b, location_t e) : beg{b}, end{e} { }
     range_t(int tid, int b, int e) : beg{utility::make_location(tid,b)}, end{utility::make_location(tid,e)} { }
     
     location_t beg;
     location_t end;
};

using ranges_t = std::deque<range_t>;
using ranges_queue_t = std::queue<range_t,ranges_t>;

// Parse a white-spaced separated list of regions
// Formats of a region range is contig:from-to
//   contig: name of contig/fragment/sequence
//   from: beginning of range, 1-based inclusive
//   to: end of range, 1-based inclusive
// Format list:
//   ctg
//   ctg:from-to
//   ctg:position
//   ctg:from-
//   ctg:-to
//   ctg:from+length

boost::optional<contig_fragments_t> parse_contig_fragments_from_regions(const std::string &text);

boost::optional<contig_fragments_t> parse_contig_fragments_from_bed(const std::string &text);

ranges_t convert_fragments_to_ranges(const contig_fragments_t& fragments, const ContigIndex& index, bool one_indexed=true);

// Parse dng region string into ranges
ranges_t parse_regions(const std::string &text, const ContigIndex& index);

// Parse slurped bed into ranges
ranges_t parse_bed(const std::string &text, const ContigIndex& index);

// Task helper functions
template<typename M>
inline
void set_regions(std::string region, const regions::ContigIndex& index, M *mpileup) {
    assert(mpileup != nullptr);
    // replace arg.region with the contents of a file if needed
    auto region_ext = io::at_slurp(region);
    if(!region.empty()) {
        if(region_ext == "bed") {
            mpileup->SetRegions(regions::parse_bed(region, index));
        } else {
            mpileup->SetRegions(regions::parse_regions(region, index));
        }
    }
}

template<typename M>
inline
void set_regions(std::string region, M *mpileup) {
    assert(mpileup != nullptr);
    // replace arg.region with the contents of a file if needed
    auto region_ext = io::at_slurp(region);
    if(!region.empty()) {
        if(region_ext == "bed") {
            if(auto f = regions::parse_contig_fragments_from_bed(region)) {
                mpileup->SetRegions(*f, false);
            } else {
                throw std::invalid_argument("Parsing of bed failed.");
            }
        } else {
            if(auto f = regions::parse_contig_fragments_from_regions(region)) {
                mpileup->SetRegions(*f, true);
            } else {
                throw std::invalid_argument("Parsing of regions failed.");
            }
        }
    }
}

} // namespace dng::regions
} // namespace dng

#endif // DNG_REGIONS_H
