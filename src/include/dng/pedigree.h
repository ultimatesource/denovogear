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


/*

##PEDNG v0.1
##
#Indiv        Dad         Mom         Sex   Samples
A1@founder    .           .           1     .
A2            .           .           2     A2a A2b
A3            A1          A2          1     =
A4@gamete     A1:0.1      .           1     =
A5@clone      A3          .           1     =
A6@ploidy=1@founder .     .           1     =

*/

#pragma once
#ifndef DNG_IO_PEDIGREE_H
#define DNG_IO_PEDIGREE_H

#include <vector>
#include <string>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <boost/optional.hpp>

#include <dng/utility.h>

namespace dng {

class Pedigree {
public:
    enum class Sex : char {
        Autosomal = 0, Male = 1, Female = 2,
        Unknown
    };

    static Sex parse_sex(const std::string &str);

    template<typename Range>
    static Pedigree parse_text(const Range &text);

    static Pedigree parse_table(const std::vector<std::vector<std::string>> &table);

    struct Member {
        std::string name;
        std::vector<std::string> tags;

        boost::optional<std::string> dad;
        boost::optional<double> dad_length;

        boost::optional<std::string> mom;
        boost::optional<double> mom_length;

        Sex sex;
 
        std::vector<std::string> samples;
    };

    using MemberTable = std::vector<Member>;

    Pedigree() = default;

    explicit Pedigree(MemberTable table) : table_{table} {
        for(MemberTable::size_type i=0; i < table_.size(); ++i) {
            auto ret = names_.emplace(table_[i].name, i);
            if(!ret.second) {
                throw std::invalid_argument("The name of a member of the pedigree is not unique: '" + table_[i].name + "'.");
            }
        }
    }

    std::size_t AddMember(Member member) {
        auto pos = table_.size();
        auto ret = names_.emplace(member.name, pos);
        if(!ret.second) {
            throw std::invalid_argument("The name of a member of the pedigree is not unique: '" + member.name + "'.");
        }
        table_.push_back(std::move(member));
        return pos;
    }

    std::size_t LookupMemberPosition(const std::string & child) const {
        auto it = names_.find(child);
        if(it == names_.end()) {
            return names_.size();
        }
        return it->second;
    }

    const Member* LookupMember(const std::string & child) const {
        auto it = names_.find(child);
        if(it == names_.end()) {
            return nullptr;
        }
        return &table_[it->second];
    }

    const Member& GetMember(std::size_t pos) const { return table_.at(pos); }

    std::size_t NumberOfMembers() const { return table_.size(); } 

    void Clear() {
        table_.clear();
        names_.clear();
    }

    const MemberTable& table() { return table_; }

private:
    MemberTable table_;
    std::unordered_map<std::string,MemberTable::size_type> names_;
};

inline
Pedigree::Sex Pedigree::parse_sex(const std::string &str) {
    static std::pair<std::string, Sex> keys[] = {
        {".", Sex::Unknown},
        {"0", Sex::Autosomal},
        {"1", Sex::Male},
        {"2", Sex::Female},
        {"male", Sex::Male},
        {"female", Sex::Female},
        {"autosomal", Sex::Autosomal}
    };
    return dng::utility::key_switch_tuple(str, keys, keys[0]).second;
}

template<typename Range>
Pedigree Pedigree::parse_text(const Range &text) {
    using namespace boost;
    using namespace std;
    // Construct the tokenizer
    // token are separated by one or more <space>s or <tab>s
    // <newline>s end the row
    auto tokens = utility::make_tokenizer_dropempty(text, "\t ", "\n");

    auto token_it = tokens.begin();
    if(token_it == tokens.end() || *token_it != "##PEDNG") {
        throw std::invalid_argument("Pedigree parsing failed; "
            "unknown pedigree format; missing '##PEDNG' header line.");
    }

    // Work through tokens and build each row
    size_t k = 0;
    bool in_comment = false;
    vector<vector<string>> string_table;
    string_table.reserve(64);
    for(; token_it != tokens.end(); ++token_it) {
        const auto & token = *token_it;
        if(token == "\n") {
            k = 0;
            in_comment = false;
            continue;
        }
        if(in_comment) {
            // skip rows that are comments
            continue;
        }
        if(k == 0) {
            string_table.emplace_back();
            if(token[0] == '#') {
                in_comment = true;
                continue;
            }
            string_table.back().reserve(8);
        }
        // Add token to the current row
        string_table.back().push_back(token);
        k += 1;
    }

    return parse_table(string_table);
 }


} // namespace dng

#endif //DNG_IO_PEDIGREE_H
