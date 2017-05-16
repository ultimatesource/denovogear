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

#include <dng/pool.h>
#include <dng/utility.h>
#include <dng/regions.h>
#include <dng/library.h>

#include <dng/hts/bam.h>
#include <dng/hts/extra.h>

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

private:
    File in_;
    utility::location_t next_loc_;
    list_type buffer_;
    int min_qlen_;
};
} // namespace detail

class BamPileup {
public:
    typedef detail::BamScan::pool_type pool_type;
    typedef detail::BamScan::list_type list_type;
    typedef detail::BamScan::node_type node_type;

    struct contig_t {
        std::string name;
        int length;
    };

    typedef std::vector<list_type> data_type;
    typedef void (callback_type)(const data_type &, utility::location_t);

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

    const libraries_t& libraries() const {
        return output_libraries_;
    }
    size_t num_libraries() const {
        return output_libraries_.names.size();
    }
    const utility::StringSet& read_groups(size_t index) const {
        return output_libraries_.read_groups[index];
    }

    const bam_hdr_t * header() const {
        return scanners_.front().file().header();
    }

    std::vector<contig_t> contigs() const {
        std::vector<contig_t> ret;
        for(auto && contig : hts::bam::contigs(header())) {
            ret.push_back({contig.first,contig.second});
        }
        return ret;
    }


protected:
    int Advance(data_type *data, utility::location_t *target_loc,
                utility::location_t *fast_forward_loc);

private:
    void ParseHeader(const char* text);

    template<typename It>
    void ParseHeaderTokens(It it, It it_last);


    std::vector<detail::BamScan> scanners_;

    pool_type pool_;

    int min_qlen_;

    struct bam_libraries_t : dng::libraries_t {
        std::vector<utility::StringSet> read_groups;
    };

    std::string lbtag_;

    bam_libraries_t input_libraries_;
    bam_libraries_t output_libraries_;

    utility::StringMap read_group_to_libraries_;

    DNG_UNIT_TEST(unittest_dng_io_bam);
};

template<typename InFile>
void BamPileup::AddFile(InFile&& f) {
    scanners_.emplace_back(std::forward<InFile>(f), min_qlen_);
    const char * text = scanners_.back().file().header()->text;
    ParseHeader(text);
}

template<typename CallBack>
void BamPileup::operator()(CallBack call_back) {
    using namespace std;
    using utility::make_location;
    using dng::regions::location_range_t;

    // data will hold our pileup information
    data_type data(num_libraries());
    location_t current_loc = 0;
    location_t fast_forward_loc = 0;

    // if we have no read_groups specified, use each file as a single group
    // TODO: check this
    if(data.empty()) {
        data.resize(scanners_.size());
    }

    // If the first file has parsed regions, use them.
    std::queue<location_range_t> region_queue;
    for(auto && r : scanners_.front().file().regions()) {
        // convert regions from bam format to dng format
        region_queue.emplace(r);
    }
    location_range_t current_reg = {0, utility::LOCATION_MAX};
    if(!region_queue.empty()) {
        current_reg = region_queue.front();
        region_queue.pop();
    }

    for(;; current_loc += 1) {
        int res = Advance(&data, &current_loc, &fast_forward_loc);
        if(res < 0) {
            break;
        } else if(res == 0) {
            continue;
        }
        // if we have advanced passed our current region, try to find the next one
        while( current_reg.end <= current_loc ) {
            if(region_queue.empty()) {
                break;
            }
            current_reg = region_queue.front();
            region_queue.pop();
        }
        // location does not overlap our region so skip it
        if( current_loc < current_reg.beg || current_reg.end <= current_loc ) {
            continue;
        }
        //if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
        // Execute callback function
        call_back(data, current_loc);
    }
}

template<typename R>
void BamPileup::SelectLibraries(R &range) {
    // Clear all output libraries
    output_libraries_ = {};

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
