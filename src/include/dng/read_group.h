/*
 * Copyright (c) 2014,2015 Reed A. Cartwright
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
#ifndef DNG_READ_GROUP_H
#define DNG_READ_GROUP_H

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/unique.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/container/flat_set.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#include <unordered_set>

#include <htslib/vcf.h>

namespace dng {

namespace rg {
struct id {};
struct lb {};
struct sm {};
struct nx {};
}

// @RG ID: str->index; index->lb_index;
// @RG LB: str->index; index->str;
// @RG SM: str->index; index->lb_indies;

namespace detail {
namespace mi = boost::multi_index;
using boost::multi_index_container;
using mi::indexed_by;
using mi::tag;
using mi::member;
using mi::identity;
using mi::ordered_non_unique;
using mi::ordered_unique;
using mi::random_access;

struct rg_t {
    std::string id;
    std::string library;
    std::string sample;
};

typedef boost::multi_index_container<rg_t, indexed_by<
ordered_unique<tag<rg::id>, member<rg_t, std::string, &rg_t::id>>,
               ordered_non_unique<tag<rg::lb>, member<rg_t, std::string, &rg_t::library>>,
               ordered_non_unique<tag<rg::sm>, member<rg_t, std::string, &rg_t::sample>>,
               random_access<tag<rg::nx>>
               >> DataBase;
}

class ReadGroups {
public:
    typedef boost::container::flat_set<std::string> StrSet;
    typedef detail::rg_t ReadGroup;
    typedef detail::DataBase DataBase;

    // Parse a list of BAM/SAM files into read groups
    template<typename InFiles>
    void ParseHeaderText(InFiles &range);

    // Parse a VCF file
    template<typename InFiles>
    void ParseSamples(InFiles &range);

    inline const StrSet &groups() const { return groups_; }
    inline const StrSet &libraries() const { return libraries_; }
    inline const StrSet &samples() const { return samples_; }
    inline const DataBase &data() const { return data_; }

    // Returns the library index based on read group sorted id
    inline std::size_t library_from_id(std::size_t n) const {
        return library_from_id_[n];
    }

    // Returns the library index based on read group order of discovery
    inline std::size_t library_from_index(std::size_t n) const {
        return library_from_index_[n];
    }

    template<typename Range>
    inline void EraseLibraries(Range &range) {
        for(auto && lib : range) {
            data_.get<rg::lb>().erase(lib);
        }
        ReloadData();
    }

protected:
    DataBase data_;

    StrSet groups_;
    StrSet libraries_;
    StrSet samples_;

    std::vector<std::size_t> library_from_id_;
    std::vector<std::size_t> library_from_index_;

    inline void ReloadData();
};

namespace rg {
inline std::size_t index(const ReadGroups::StrSet &set,
                         const std::string &query) {
    auto it = set.find(query);
    return (it != set.end()) ? static_cast<std::size_t>(it - set.begin()) : -1;
}
}

inline void ReadGroups::ReloadData() {
    // clear data vectors
    groups_.clear();
    libraries_.clear();
    samples_.clear();

    // construct data vectors
    groups_.reserve(data_.size());
    libraries_.reserve(data_.size());
    samples_.reserve(data_.size());
    for(auto && a : data_) {
        groups_.insert(a.id);
        libraries_.insert(a.library);
        samples_.insert(a.sample);
    }

    // remember the order the libraries were discovered
    library_from_index_.reserve(groups_.size());
    for(std::size_t u = 0; u < data_.get<rg::nx>().size(); ++u) {
        library_from_index_.push_back(rg::index(libraries_,
                                                data_.get<rg::nx>()[u].library));
    }

    // connect groups to libraries and samples
    library_from_id_.resize(groups_.size());
    for(auto && a : data_) {
        library_from_id_[rg::index(groups_, a.id)] = rg::index(libraries_, a.library);

    }
}

template<typename InFiles>
void ReadGroups::ParseHeaderText(InFiles &range) {
    // iterate through each file in the range
    for(auto && f : range) {
        // get header text and continue on failure
        const char *text = f.header()->text;
        if(text == nullptr) {
            continue;
        }

        // enumerate over read groups
        for(text = strstr(text, "@RG\t"); text != nullptr;
                text = strstr(text, "@RG\t")) {
            text += 4; // skip @RG\t
            // parse the @RG line as tab-separated key:value pairs
            // use a map to separate the parsing logic from the rg_t struct
            const char *k = text, *v = k, *p = k;
            std::map<std::string, std::string> tags;
            for(; *p != '\n' && *p != '\0'; ++p) {
                if(*p == ':') {
                    v = p + 1;    // make v point the the beginning of the value
                }
                // v-1 will point to the end of the key
                else if(*p == '\t') {
                    if(k < v) {
                        // if we found a key:value pair, add it to the map
                        // first one wins
                        tags.emplace(std::string(k, v - 1), std::string(v, p));
                    }
                    // move k and v to the next character
                    k = v = p + 1;
                }
            }
            // make sure the insert the last item if needed
            if(k < v) {
                tags.emplace(std::string(k, v - 1), std::string(v, p));
            }

            // continue to next @RG if this @RG does not have an ID
            auto it = tags.find("ID");
            if(it == tags.end() || it->second.empty()) {
                continue;
            }
            // check to see if this is a duplicate
            // TODO: Raise warning/error if this is a collision of different sm/lb data
            if(data_.find(it->second) != data_.end()) {
                continue;
            }
            ReadGroup val{it->second};

            // sample tag
            if((it = tags.find("SM")) != tags.end()) {
                val.sample = std::move(it->second);
            } else {
                val.sample = val.id;
            }
            // library tag
            it = tags.find("LB");
            val.library = val.sample;
            if(it != tags.end()) {
                // Prevent collision of LB tags from different samples by prepending sample tag
                if(val.library != it->second) {
                    val.library += '-' + it->second;
                }
            } else if(val.library != val.id) {
                val.library += '-' + val.id;
            }
            data_.insert(std::move(val));
        }
    }
    ReloadData();
}

/**
 * Parses an VCF file into read groups, sample names and libraries
 * fname - file path and name string for input vcf file
 */
template<typename InFile>
void ReadGroups::ParseSamples(InFile &file) {
    // Read through samples from the header
    auto samples = file.samples();
    for(int a = 0; a < samples.second; ++a) {
        // Each sample column in a vcf file typically corresponds to an @RG SM: tag.
        // However, we prefer to work with @RG LB: tags as our base data.
        // In order to extract library information from VCF data, we support column
        // ids in the following format SAMPLE:LIBRARY.

        const char *sm = samples.first[a];
        ReadGroup val{sm, sm};

        const char *p = strchr(sm, ':');
        if(p == nullptr || *(p + 1) == '\0') {
            val.sample = sm;
        } else {
            val.sample.assign(sm, p);
        }
        data_.insert(std::move(val));
    }
    ReloadData();
}

} // namespace dng

#endif
