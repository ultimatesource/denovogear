/*
 * Copyright (c) 2014-2016 Reed A. Cartwright
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

#include <dng/read_group.h>

using namespace std;

namespace dng {

void ReadGroups::ReloadData() {
    // clear data vectors
    groups_.clear();
    libraries_.clear();
    samples_.clear();
    library_from_index_.clear();
    library_from_id_.clear();

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

void ReadGroups::ParseHeaderText(const char* text, const std::string &lbtag) {
    // Get tag to identify libraries, LB by default.
    const std::string &lbtag_str = (lbtag.empty() ? "LB" : lbtag);

    if(text == nullptr) {
        return;
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
            val.sample = it->second;//std::move(it->second);
        } else {
            val.sample = val.id;
        }

        // library tag
        it = tags.find(lbtag_str);
        if(it == tags.end()) {
            // Will throw an error if any @RG header line is missing tag.
            throw std::runtime_error("An @RG header is missing " + lbtag_str + " tag." +
                                     "  Please fix or use option \'--lbtag\' to use another tag to identify libraries.");
        } else {
            // Unless using SM tag to group different libraries, prepend SM tag for readability
            if(lbtag_str == "SM" || val.sample == it->second) {
                val.library = it->second;
            } else {
                val.library = val.sample + "-" + it->second;
            }
        }

        data_.insert(std::move(val));
    }
    ReloadData();
}

}; // namespace dng