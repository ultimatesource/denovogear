/*
 * Copyright (c) 2014-2017 Reed A. Cartwright
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
#ifndef DNG_IO_BAM_H
#define DNG_IO_BAM_H

#include <vector>
#include <queue>
#include <unordered_map>
#include <utility>

#include <boost/optional.hpp>

#include <dng/pool.h>
#include <dng/utility.h>
#include <dng/regions.h>
#include <dng/library.h>

#include <dng/hts/bam.h>

#include <dng/detail/unit_test.h>

namespace dng {
namespace io {

struct bam_record_t : public PoolNode {
    hts::bam::Alignment aln; // sequence record
    utility::location_t beg; // 0-based left-most edge, [beg,end)
    utility::location_t end; // 0-based right-most edge, [beg,end)
    utility::location_t pos; // current position of the pileup in this query/read
    bool is_missing; // is there no base call at the pileup position?

    hts::bam::cigar_t cigar; // cigar
    hts::bam::data_t seq;    // encoded sequence
    hts::bam::data_t qual;   // quality scores

   inline uint8_t base() const {
        assert(0 <= pos && pos < (qual.second-qual.first));
        return bam_seqi(seq.first, pos);
    }

   inline uint8_t base_qual() const {
        assert(0 <= pos && pos < (qual.second-qual.first));
        return qual.first[pos];
    }
};

namespace detail {
using BamPool = IntrusivePool<bam_record_t>; 

class BamScan {
public:
    typedef hts::bam::File File;
    typedef detail::BamPool pool_type;
    typedef pool_type::list_type list_type;
    typedef pool_type::node_type node_type;

    explicit BamScan(File in, int min_qlen = 0) : in_(std::move(in)), next_loc_{0},
        min_qlen_{min_qlen} {}

    list_type operator()(utility::location_t target_loc, pool_type &pool);

    location_t next_loc() const { return next_loc_; }

    const File& file() const { return in_; }

    void SetRegion(const regions::range_t &region) {
        int tid = utility::location_to_contig(region.beg);
        assert( tid == utility::location_to_contig(region.end));
        int beg = utility::location_to_position(region.beg);
        int end = utility::location_to_position(region.end);
        in_.SetRegion(tid, beg, end);        
    }

private:
    File in_;
    utility::location_t next_loc_;
    list_type buffer_;
    int min_qlen_;
};
} // namespace detail

class BamPileup {
public:
    using pool_type = detail::BamScan::pool_type;
    using list_type = detail::BamScan::list_type;
    using node_type = detail::BamScan::node_type;

    using data_type = std::vector<list_type>;
    using callback_type = void(const data_type &, utility::location_t);

    template<typename CallBack>
    void operator()(CallBack func);

    BamPileup(int min_qlen = 0, std::string lbtag = "LB") : pool_{4048}, min_qlen_{min_qlen},
        lbtag_{std::move(lbtag)} {
        if(lbtag_.empty()) {
            lbtag_ = "LB";
        }
    }

    template<typename InFile>
    void AddFile(InFile&& f);

    template<typename R>
    void SelectLibraries(R &range);

    void ResetLibraries();

    void SetRegions(regions::ranges_t regions) {
        regions_ = regions::ranges_queue_t(std::move(regions));
    }

    const libraries_t& libraries() const {
        return output_libraries_;
    }
    size_t num_libraries() const {
        return output_libraries_.names.size();
    }
    const utility::StringSet& read_groups(size_t index) const {
        return output_libraries_.read_groups[index];
    }

    // retrieves the header from the first file
    const bam_hdr_t * header() const {
        return scanners_.front().file().header();
    }

    // retrieves the contigs from the first file
    std::vector<regions::contig_t> contigs() const {
        std::vector<regions::contig_t> ret;
        for(auto && contig : hts::bam::contigs(header())) {
            ret.push_back({std::string{contig.first}, contig.second});
        }
        return ret;
    }


private:
    int Advance(utility::location_t *target_location);

    // Data Used for Pileup
    data_type data_; // Store pileup
    utility::location_t fast_forward_location_; // Next location off the scanners

    void ParseHeader(const char* text);

    template<typename It>
    void ParseHeaderTokens(It it, It it_last);

    boost::optional<regions::range_t> LoadNextRegion();

    std::vector<detail::BamScan> scanners_;

    regions::ranges_queue_t regions_;

    pool_type pool_;

    int min_qlen_;

    struct bam_libraries_t : dng::libraries_t {
        std::vector<utility::StringSet> read_groups;
    };

    std::string lbtag_;

    bam_libraries_t input_libraries_;
    bam_libraries_t output_libraries_;

    utility::StringMap read_group_to_libraries_;

    DNG_UNIT_TEST_CLASS(unittest_dng_io_bam);
};

template<typename InFile>
void BamPileup::AddFile(InFile&& f) {
    scanners_.emplace_back(std::forward<InFile>(f), min_qlen_);
    const char * text = scanners_.back().file().header()->text;
    ParseHeader(text);
}

inline
boost::optional<regions::range_t> BamPileup::LoadNextRegion() {
    using regions::range_t;

    if(regions_.empty()) {
        return boost::none;
    }
    range_t region = regions_.front();
    regions_.pop();
    for(auto &&s : scanners_) {
        s.SetRegion(region);
    }

    return region;
}

template<typename CallBack>
void BamPileup::operator()(CallBack call_back) {
    using namespace std;
    using utility::make_location;
    using dng::regions::range_t;

    // data will hold our pileup information
    data_.clear();
    data_.resize(num_libraries());
    fast_forward_location_ = 0;

    range_t current_reg = {0, utility::LOCATION_MAX};
    if(auto reg = LoadNextRegion()) {
        current_reg = *reg;
    }
    utility::location_t current_location = current_reg.beg;

    for(;; current_location += 1) {
        int res = Advance(&current_location);
        while( res < 0 ) {
            // We are out of reads, try loading the next region
            if(auto reg = LoadNextRegion()) {
                current_reg = *reg;
                current_location = current_reg.beg;
                res = Advance(&current_location);
            } else {
                // we are out of regions
                break;
            }
        }
        if( res < 0 ) {
            // we are out of data
            break;
        } else if(res == 0) {
            // We are not out of reads, but nothing overlaps current_loc
            continue;            
        }

        // location does not overlap our region so skip it
        if( current_location < current_reg.beg || current_reg.end <= current_location ) {
            continue;
        }
        // Execute callback function
        call_back(data_, current_location);
    }
}

template<typename R>
void BamPileup::SelectLibraries(R &range) {
    // Clear all output libraries
    output_libraries_ = {};
    read_group_to_libraries_.clear();

    // For every library in range, try to find it in input_libraries_
    size_t k=0;
    for(auto it = boost::begin(range); it != boost::end(range); ++it) {
        auto pos = utility::find_position(input_libraries_.names, *it);
        if(pos == input_libraries_.names.size()) {
            // Do nothing if library was not found.
            continue;
        }
        output_libraries_.names.push_back(input_libraries_.names[pos]);
        output_libraries_.samples.push_back(input_libraries_.samples[pos]);
        output_libraries_.read_groups.push_back(input_libraries_.read_groups[pos]);
        for(auto && a : input_libraries_.read_groups[pos]) {
            read_group_to_libraries_.emplace(a,k);
        }
        ++k;
    }

}

} //namespace io
} //namespace dng

#endif //DNG_IO_BAM_H
