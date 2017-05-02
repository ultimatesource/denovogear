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

#include <dng/io/bam.h>
#include <dng/cigar.h>
#include <dng/utility.h>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/metafunctions.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/sort.hpp>

#include <boost/algorithm/string/predicate.hpp>

using namespace dng;
using namespace dng::io;
using namespace dng::io::detail;

// Advance starts a pileup procedure at pos.
// If there are no reads at pos, it forwards to the next pos and updates
// current_pos
int BamPileup::Advance(data_type *data, utility::location_t *target_loc,
                       utility::location_t *fast_forward_loc) {
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
            location_t q = cigar::target_to_query(*target_loc, it->beg, it->cigar);
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
        location_t next_loc = utility::LOCATION_MAX;
        std::size_t k = 0;
        // Iterate over input files
        for(auto &scanner : scanners_) {
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
                if(!read_group_to_libraries_.empty()) {
                    auto it = read_group_to_libraries_.find(reinterpret_cast<const char *>(rg + 1));
                    if(it == read_group_to_libraries_.end()) {
                        pool_.Free(p); // drop unknown RG's
                        continue;
                    }
                    index = it->second;
                }
                // process cigar string
                location_t q = cigar::target_to_query(*target_loc, p->beg, p->cigar);
                p->pos = cigar::query_pos(q);
                p->is_missing = cigar::query_del(q);

                // push read onto correct RG
                (*data)[index].push_back(*p);
                fast_forward = false;
            }
            ++k;
        }
        *fast_forward_loc = next_loc;
    }
    // fast_forward will be true if no reads overlap our current position
    // return -1 if we are out of reads and at the end of our reference range
    if(fast_forward && *fast_forward_loc == utility::LOCATION_MAX) {
        return -1;
    }
    // return 1 if we have read data
    // return 0 if we don't have read data at this location but we may have read data in the future 
    return (fast_forward ? 0 : 1);
}

BamScan::list_type
BamScan::operator()(utility::location_t target_loc, pool_type &pool) {
    using utility::make_location;
    // Reads from in_ until it encounters the first read
    // that is right of pos.
    // TODO: check for proper read ordering???
    while(next_loc_ <= target_loc) {
        node_type *p = pool.Malloc();
        do {
            // Try to grab a read, if not return current buffer
            if(in_(&p->aln) < 0) {
                pool.Free(p);
                next_loc_ = std::numeric_limits<location_t>::max();
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

void BamPileup::ParseHeader(const char* text) {
    if(text == nullptr) {
        return;
    }
    // Tokenize header
    auto tokens = utility::make_tokenizer_dropempty(text);
    // Parse Tokens
    ParseHeaderTokens(tokens.begin(), tokens.end());    
}

template<typename It>
void BamPileup::ParseHeaderTokens(It it, It it_last) {
    using boost::algorithm::starts_with;

    bool needs_updating = false;

    std::string prefix = lbtag_ + ':';
    for(; it != it_last; ++it) {
        if(*it == "@RG") {
            std::string name, id, sample;
            for(++it; it != it_last && *it != "\n"; ++it) {
                if(starts_with(*it, "ID:")) {
                    id = it->substr(3);
                } else if(starts_with(*it, prefix)) {
                    name = it->substr(prefix.size());
                } else if(starts_with(*it, "SM:")) {
                    sample = it->substr(3);
                }
            }
            if(id.empty()) {
                throw std::runtime_error("An @RG header line is missing an 'ID' tag.");
            }
            if(name.empty()) {
                name = id;
            }
            if(sample.empty()) {
                sample = name;
            }
            auto pos = utility::find_position(input_libraries_.names, name);
            if(pos == input_libraries_.names.size()) {
                input_libraries_.names.push_back(std::move(name));
                input_libraries_.samples.push_back(std::move(sample));
                input_libraries_.read_groups.push_back(utility::StringSet{{id}});
                needs_updating = true;
            } else {
                if(input_libraries_.samples[pos] != sample) {
                    throw std::runtime_error("Multiple sample names are defined for read groups with library " + prefix + name + ".");
                }
                needs_updating |= input_libraries_.read_groups[pos].insert(id).second;
            }
        }
    }
    if(needs_updating) {
        ResetLibraries();
    }
}


void BamPileup::ResetLibraries() {
    output_libraries_ = input_libraries_;
    
    for(size_t i = 0; i < output_libraries_.read_groups.size(); ++i) {
        for(auto && a : output_libraries_.read_groups[i]) {
            read_group_to_libraries_.emplace(a,i);
        }
    }
}
