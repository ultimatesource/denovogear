/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
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
#include <unordered_map>

#include <dng/hts/bam.h>
#include <htslib/faidx.h>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/metafunctions.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/sort.hpp>

#include <dng/pool.h>
#include <dng/fileio.h>
#include <dng/cigar.h>
#include <dng/read_group.h>

namespace dng {

inline uint64_t make_location(int t, int p) {
    return (static_cast<uint64_t>(t) << 32) | p;
}
inline int location_to_target(uint64_t u) {
    return static_cast<int>(u >> 32);
}
inline int location_to_position(uint64_t u) {
    return static_cast<int>(u & 0x7FFFFFFF);
}

namespace pileup {

struct Node : public PoolNode {
    hts::bam::Alignment aln;   // sequence record
    uint64_t beg; // 0-based left-most edge, [beg,end)
    uint64_t end; // 0-based right-most edge, [beg,end)
    uint64_t pos; // current position of the pileup in this query/read
    bool is_missing; // is there no base call at the pileup position?

    hts::bam::cigar_t cigar; // cigar
    hts::bam::data_t seq;    // encoded sequence
    hts::bam::data_t qual;   // quality scores
};

namespace detail {
typedef IntrusivePool<Node> NodePool;
}

typedef detail::NodePool::list_type NodeList;

namespace detail {

template<typename InFile>
class BamScan {
public:
    typedef NodeList list_type;
    typedef Node node_type;
    typedef detail::NodePool pool_type;

    explicit BamScan(InFile &in, int min_qlen = 0) : in_{in}, next_loc_{0},
        min_qlen_{min_qlen} {}

    list_type operator()(uint64_t target_loc, pool_type &pool) {
        // Reads from in_ until it encounters the first read
        // that is right of pos.
        // TODO: check for proper read ordering???
        while(next_loc_ <= target_loc) {
            node_type *p = pool.Malloc();
            do {
                // Try to grab a read, if not return current buffer
                if(in_(&p->aln) < 0) {
                    pool.Free(p);
                    next_loc_ = std::numeric_limits<uint64_t>::max();
                    return std::move(buffer_);
                }
                p->cigar = p->aln.cigar();

                if(min_qlen_ > 0 && cigar::query_length(p->cigar) < min_qlen_) {
                    continue;
                }

                p->beg = make_location(p->aln.target_id(), p->aln.position());
                // update right-most position in the read
                p->end = p->beg + cigar::target_length(p->cigar);
            } while(p->end <= target_loc);

            // cache pointers to data elements
            p->seq = p->aln.seq();
            p->qual = p->aln.seq_qual();
            // TODO: Adjust for Illumina 1.3 quality

            // update the left-most position of the most recently read read.
            next_loc_ = p->beg;
            // save read
            buffer_.push_back(*p);
        }
        // Return all but the last read.
        if(buffer_.size() <= 1) {
            return list_type{};
        }
        list_type ret;
        ret.splice(ret.end(), buffer_, buffer_.begin(), --buffer_.end());
        return ret;
    }

    uint64_t next_loc() const { return next_loc_; }

private:
    uint64_t next_loc_;
    InFile &in_;
    list_type buffer_;
    int min_qlen_;
};

} // namespace detail

class BamPileup {
public:
    typedef NodeList list_type;
    typedef Node node_type;
    typedef detail::NodePool pool_type;

    typedef std::vector<list_type> data_type;
    typedef void (callback_type)(const data_type &, uint64_t);

    template<typename InFiles, typename Func>
    void operator()(InFiles &range, Func func);

    template<typename RG>
    BamPileup(const RG &rg, int min_qlen = 0) : pool_{4048}, read_groups_{rg},
        min_qlen_{min_qlen} {

    }

    BamPileup() { }

protected:
    template<typename Scanners>
    int Advance(Scanners &range, data_type *data, uint64_t *target_loc,
                uint64_t *fast_forward_loc);

    template<typename STR>
    std::size_t ReadGroupIndex(const STR &s) {
        auto it = read_groups_.find(s);
        return (it == read_groups_.end()) ? -1
               : static_cast<std::size_t>(it - read_groups_.begin());
    }

private:
    pool_type pool_;

    ReadGroups::StrSet read_groups_;

    int min_qlen_;
};

template<typename InFiles, typename Func>
void BamPileup::operator()(InFiles &range, Func func) {
    using namespace std;
    using namespace fileio;

    // encapsulate input files into scanners
    vector<detail::BamScan<typename boost::range_value<InFiles>::type>> scanners;
    for(auto it =  boost::begin(range); it != boost::end(range); ++it) {
        scanners.emplace_back(*it, min_qlen_);
    }

    // type erase our callback function
    function<callback_type> call_back(func);

    // data will hold our pileup information
    data_type data(read_groups_.size());
    uint64_t current_loc = 0;
    uint64_t fast_forward_loc = 0;

    // if we have no read_groups specified, use each file as a single group
    // TODO: check this
    if(data.empty()) {
        data.resize(scanners.size());
    }

    // If there is a parsed region, use it.
    // TODO: check to see if the regions are all the same???
    uint64_t beg_loc = 0, end_loc = std::numeric_limits<uint64_t>::max();
    if(boost::begin(range)->iter() != nullptr) {
        beg_loc = make_location(boost::begin(range)->iter()->tid,
                                boost::begin(range)->iter()->beg);
        end_loc = make_location(boost::begin(range)->iter()->tid,
                                boost::begin(range)->iter()->end);
    }
    for(;; current_loc += 1) {
        int res = Advance(scanners, &data, &current_loc, &fast_forward_loc);
        if(res < 0) {
            break;
        }
        if(res == 0 || beg_loc > current_loc || current_loc >= end_loc) {
            continue;
        }
        //if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
        // Execute callback function
        call_back(data, current_loc);
    }
}

// Advance starts a pileup procedure at pos.
// If there are no reads at pos, it forwards to the next pos and updates
// current_pos

template<typename Scanners>
int BamPileup::Advance(Scanners &range, data_type *data, uint64_t *target_loc,
                       uint64_t *fast_forward_loc) {
    using namespace std;
    // enumerate through existing data set, updating location of reads and
    // purge reads that have expired
    bool fast_forward = true;
    for(auto &d : *data) {
        auto it = d.begin();
        while(it != d.end()) {
            if(it->end <= *target_loc) {
                node_type *p = &(*it);
                d.erase(it++);
                pool_.Free(p);
                continue;
            }
            uint64_t q = cigar::target_to_query(*target_loc, it->beg, it->cigar);
            it->pos = cigar::query_pos(q);
            it->is_missing = cigar::query_del(q);
            ++it;
        }
        // we want to fast_forward if there is nothing to output
        fast_forward = (fast_forward && d.empty());
    }
    if(fast_forward) {
        *target_loc = *fast_forward_loc;
    }
    if(*target_loc >= *fast_forward_loc) {
        uint64_t next_loc = numeric_limits<uint64_t>::max();
        std::size_t k = 0;
        for(auto it = boost::begin(range); it != boost::end(range); ++it, ++k) {
            auto &scanner = *it;
            // Scan reads from file until target_loc
            list_type new_reads = scanner(*target_loc, pool_);
            // Update the minimum position of the next read
            next_loc = std::min(next_loc, scanner.next_loc());
            // process read_groups
            while(!new_reads.empty()) {
                node_type *p = &new_reads.front();
                new_reads.pop_front();
                uint8_t *rg = p->aln.aux_get("RG");
                std::size_t index = k;
                if(!read_groups_.empty()) {
                    index = ReadGroupIndex(reinterpret_cast<const char *>(rg + 1));
                    if(index == -1) {
                        pool_.Free(p); // drop unknown RG's
                        continue;
                    }
                }
                // process cigar string
                uint64_t q = cigar::target_to_query(*target_loc, p->beg, p->cigar);
                p->pos = cigar::query_pos(q);
                p->is_missing = cigar::query_del(q);

                // push read onto correct RG
                (*data)[index].push_back(*p);
                fast_forward = false;
            }
        }
        *fast_forward_loc = next_loc;
    }
    if(fast_forward && *fast_forward_loc == numeric_limits<uint64_t>::max()) {
        return -1;
    }
    return (fast_forward ? 0 : 1);
}
} //namespace pileup

using pileup::BamPileup;
} //namespace dng

#endif //DNG_PILEUP_H
