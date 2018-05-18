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

#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <boost/range/algorithm_ext/iota.hpp>

#include <dng/pool.h>
#include <dng/utility.h>
#include <dng/regions.h>
#include <dng/library.h>
#include <dng/depths.h>
#include <dng/utility.h>
#include <dng/seq.h>

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

    BamScan(BamScan&&) = default;

    list_type operator()(utility::location_t target_loc, pool_type &pool);

    location_t next_loc() const { return next_loc_; }

    const File& file() const { return in_; }

    void SetRegion(const regions::range_t &region) {
        int tid = utility::location_to_contig(region.beg);
        assert( tid == utility::location_to_contig(region.end));
        int beg = utility::location_to_position(region.beg);
        int end = utility::location_to_position(region.end);
        in_.SetRegion(tid, beg, end);
        next_loc_ = 0;
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

    struct Alleles; // functor class for calculating depths from data_type

    template<typename CallBack>
    void operator()(CallBack func);

    BamPileup(int min_qlen = 0, std::string lbtag = "LB") : pool_{4048}, min_qlen_{min_qlen},
        lbtag_{std::move(lbtag)} {
        if(lbtag_.empty()) {
            lbtag_ = "LB";
        }
    }
    BamPileup(BamPileup&&) = default;
    BamPileup& operator=(BamPileup&&) = default;

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

    template<typename A>
    static BamPileup open_and_setup(const A& arg);

private:
    int Advance(regions::range_t *target_range);

    void ClearData() {
        for(auto & d : data_) {
            auto it = d.begin();
            while(it != d.end()) {
                node_type *p = &(*it);
                d.erase(it++);
                pool_.Free(p);
            }
        }    
    }

    boost::optional<regions::range_t> LoadNextRegion();

    template<typename It>
    void ParseHeaderTokens(It it, It it_last);

    void ParseHeader(const char* text);


    struct bam_libraries_t : dng::libraries_t {
        std::vector<utility::StringSet> read_groups;
    };

    // Data Used for Pileup
    data_type data_; // Store pileup
    utility::location_t next_scanner_location_; // Next location off the scanners

    std::vector<detail::BamScan> scanners_;

    regions::ranges_queue_t regions_;

    pool_type pool_;

    int min_qlen_;

    std::string lbtag_;

    bam_libraries_t input_libraries_;
    bam_libraries_t output_libraries_;

    utility::StringMap read_group_to_libraries_;

    DNG_UNIT_TEST_CLASS(unittest_dng_io_bam);
};

template<typename A>
BamPileup BamPileup::open_and_setup(const A& arg) {

    BamPileup mpileup{arg.min_qlen, arg.rgtag};
    
    for(auto && str : arg.input) {
        auto file = utility::extract_file_type(str);

        hts::bam::File input{file.path.c_str(), "r", arg.fasta.c_str(), arg.min_mapqual, arg.header.c_str()};
        if(!input.is_open()) {
            throw std::runtime_error("Unable to open bam/sam/cram input file '" + str + "' for reading.");
        }
        mpileup.AddFile(std::move(input));
    }

    // Load contigs into an index
    regions::ContigIndex index;
    for(auto && a : mpileup.contigs()) {
        index.AddContig(std::move(a));
    }
    regions::set_regions(arg.region, index, &mpileup);

    return mpileup;
}

template<typename InFile>
void BamPileup::AddFile(InFile&& f) {
    scanners_.emplace_back(std::forward<InFile>(f), min_qlen_);
    const char * text = scanners_.back().file().header()->text;
    ParseHeader(text);
}

inline
boost::optional<regions::range_t> BamPileup::LoadNextRegion() {
    using regions::range_t;

    next_scanner_location_ = 0;
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

    // Resize data and free any existing nodes
    data_.resize(num_libraries());
    ClearData();

    range_t current_range = {0, utility::LOCATION_MAX};
    if(auto reg = LoadNextRegion()) {
        current_range = *reg;
    }
    for(;;) {
        Advance(&current_range);
        assert(current_range.beg <= current_range.end);
        while(current_range.beg == current_range.end) {
            // We are out of reads, try loading the next region
            ClearData();
            if(auto reg = LoadNextRegion()) {
                current_range = *reg;
                Advance(&current_range);
            } else {
                // we are out of regions
                break;
            }
        }
        if(current_range.beg == current_range.end) {
            // we are out of data
            break;
        }
        // Execute callback function
        call_back(data_, current_range.beg);
        current_range.beg += 1;
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

// This functor will convert aligned reads to read_depths
struct BamPileup::Alleles {
    using read_depths_t = dng::pileup::allele_depths_t;
    using data_type = BamPileup::data_type;
    using filter_signature = bool(const data_type::value_type &);

    Alleles(size_t num_libraries);

    template<typename F = filter_signature>
    const read_depths_t& operator()(const data_type &data, size_t ref_index, F filter
        = [](const data_type::value_type &) { return true; } );

    const std::string& alleles_str() const {
        buffer.clear();
        size_t sz = sorted.shape()[1];
        for(size_t u=0;u<sz;++u) {
            if(u > 0) {
                buffer += ',';
            }
            buffer += seq::indexed_char(indexes[u]);
        }
        return buffer;
    }

    // temporary data
    read_depths_t unsorted;
    read_depths_t sorted;
    std::vector<int> indexes;

private:
    mutable std::string buffer;
};

// Allocate workspace based on number of libraries
inline BamPileup::Alleles::Alleles(size_t num_libraries) :
    unsorted{utility::make_array(num_libraries,5u)},
    sorted{utility::make_array(num_libraries,5u)},
    indexes(5)
{
    /*noop*/;
}

template<typename F>
inline
const BamPileup::Alleles::read_depths_t&
BamPileup::Alleles::operator()(const data_type &data, size_t ref_index, F filter) {
    // reset all depth counters
    std::array<int,5> total_unsorted{0,0,0,0,0};
    std::fill_n(unsorted.data(), unsorted.num_elements(), 0);

    // pileup on read counts
    for(std::size_t u = 0; u < data.size(); ++u) {
        for(auto && r : data[u]) {
            if(filter(r)) {
                continue;
            }
            std::size_t base = seq::base_index(r.base());
            assert(unsorted[u][base] < 65535);
            unsorted[u][base] += 1;
            total_unsorted[base] += 1;
        }
    }

    // sort read counts
    indexes.resize(total_unsorted.size());
    boost::iota(indexes,0);
    boost::sort(indexes, [&total_unsorted,ref_index](int l, int r) {
        return l == ref_index || total_unsorted[l] > total_unsorted[r];
    });
    size_t sz = indexes.size();
    for(; sz > 0 && total_unsorted[indexes[sz-1]] == 0; --sz) {
        /*noop*/;
    }
    sorted.resize(utility::make_array(sorted.size(), sz));
    for(size_t i=0;i<sorted.size();++i) {
        for(size_t u=0;u<sz;++u) {
            sorted[i][u] = unsorted[i][indexes[u]];
        }
    }
    return sorted;
}

} //namespace io
} //namespace dng

#endif //DNG_IO_BAM_H
