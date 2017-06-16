/*
 * Copyright (c) 2017 Reed A. Cartwright
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
#ifndef DNG_IO_PEDIGREE_H
#define DNG_IO_PEDIGREE_H

#include <vector>
#include <string>
#include <array>
#include <unordered_map>
#include <unordered_set>

#include <dng/utility.h>

namespace dng {

class Pedigree {
public:
    enum class Sex : char {
        Unknown = 0, Male = 1, Female = 2
    };

    static Sex parse_sex(std::string &str) {
        static std::pair<std::string, Sex> keys[] = {
            {"0", Sex::Unknown},
            {"1", Sex::Male},
            {"2", Sex::Female},
            {"male", Sex::Male},
            {"female", Sex::Female},
            {"unknown", Sex::Unknown}
        };
        return dng::utility::key_switch_tuple(str, keys, keys[0]).second;
    }

    template<typename Range>
    static Pedigree parse_text(const Range &text);

    struct Member {
        std::string child;
        std::string dad;
        std::string mom;
        Sex sex;
        std::string attribute;
    };

    typedef std::vector<Member> MemberTable;

    Pedigree() = default;

    explicit Pedigree(MemberTable table) : table_{table} {
        for(MemberTable::size_type i=0; i < table_.size(); ++i) {
            auto ret = names_.emplace(table_[i].child, i);
            if(!ret.second) {
                throw std::invalid_argument("The name of a member of the pedigree is not unique: '" + table_[i].child + "'.");
            }
        }
    }

    std::size_t AddMember(Member member) {
        auto pos = table_.size();
        auto ret = names_.emplace(member.child, pos);
        if(!ret.second) {
            throw std::invalid_argument("The name of a member of the pedigree is not unique: '" + member.child + "'.");
        }
        table_.push_back(std::move(member));
        return pos;
    }

    std::size_t AddMember(std::string child, std::string dad, std::string mom,
        Sex sex, std::string attribute) {
        return AddMember({std::move(child), std::move(dad), std::move(mom), sex, std::move(attribute)});
    }

    std::size_t LookupMemberPosition(const std::string & child) const {
        auto it = names_.find(child);
        if(it == names_.end()) {
            return names_.size();
        }
        return it->second;
    }

    MemberTable::const_iterator LookupMember(const std::string & child) const {
        auto it = names_.find(child);
        if(it == names_.end()) {
            return table_.end();
        }
        return table_.begin()+it->second;
    }

    const Member& GetMember(std::size_t pos) const { return table_.at(pos); }

    std::size_t NumberOfMembers() const { return table_.size(); } 

    void Clear() {
        table_.clear();
        names_.clear();
    }

private:
    MemberTable table_;
    std::unordered_map<std::string,MemberTable::size_type> names_;
};

template<typename Range>
Pedigree Pedigree::parse_text(const Range &text) {
    using namespace boost;
    using namespace std;
    // Construct the tokenizer
    auto tokens = utility::make_tokenizer(text);

    // Work through tokens and build each row
    std::size_t k = 0;
    vector<std::array<string, 6>> string_table;
    string_table.reserve(64);
    for(auto tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
        if(*tok_iter == "\n") {
            k = 0;
            continue;
        } else if(k >= 6) {
            // ignore columns above 6
            continue;
        } else if(k == 0) {
            string_table.emplace_back();
        }
        // Add token to the current row
        string_table.back()[k] = *tok_iter;
        k += 1;
    }

    // check for unique child names
    unordered_set<string> names;
    bool use_family = false;
    for(auto &&a : string_table) {
        if(!a[1].empty() && !names.insert(a[1]).second) {
            use_family = true;
            break;
        }
    }
    // Build Pedigree Object
    Pedigree ret;
    for(auto &&a : string_table) {
        if(a[1].empty()) {
            continue;
        }
        if(use_family) {
            string child = a[0] + '_' + a[1];
            string dad = a[0] + '_' + a[2];
            string mom = a[0] + '_' + a[3];
            ret.AddMember(std::move(child), std::move(dad), std::move(mom),
                Pedigree::parse_sex(a[4]), std::move(a[5]));
        } else {
            ret.AddMember(std::move(a[1]), std::move(a[2]), std::move(a[3]),
                Pedigree::parse_sex(a[4]), std::move(a[5]));
        }
    }
    return ret;
}


} // namespace dng

#endif //DNG_IO_PEDIGREE_H
