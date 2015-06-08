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
#ifndef DNG_IO_PED_H
#define DNG_IO_PED_H

#include <iostream>
#include <vector>
#include <array>
#include <deque>
#include <string>
#include <map>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

#include <boost/tokenizer.hpp>

#include <dng/utilities.h>

namespace dng {
namespace io {

class Pedigree {
public:
    typedef boost::multi_index_container<std::string,
            boost::multi_index::indexed_by<
            boost::multi_index::random_access<>,
            boost::multi_index::ordered_unique<
            boost::multi_index::identity<std::string>>
            >> NameContainer;

    enum class Gender : char {
        Unknown = 0, Male = 1, Female = 2
    };

    struct Member {
        std::size_t fam;
        std::size_t child;
        std::size_t dad;
        std::size_t mom;
        Gender sex;
        std::string sample_tree;
    };

    typedef std::vector<Member> MemberTable;

    // How many members, including the dummy, are in the pedigree.
    std::size_t member_count() const {
        return names_.size();
    }

    // Fetch the name of a member of the pedigree
    const std::string &name(std::size_t id) const {
        return names_[id];
    }

    // Given the name of an individual return its id.
    // If the name is not valid, return the id of the 0-th
    // individual.
    std::size_t id(const std::string &name) const {
        auto it = names_.get<1>().find(name);
        if(it == names_.get<1>().end()) {
            return std::size_t(0);
        }
        return names_.project<0>(it) - names_.begin();
    }

    const MemberTable &table() const {
        return table_;
    }

    // Parse a string-like object into a pedigree
    // Insert additional founders such that every child either is a
    //     founder or has both parents in the tree.
    // Sort the individuals starting with the founders, and then working
    //     the way down.
    // TODO: maybe we can just keep pointers to the delimiters in memory
    // TODO: add comment support
    // TODO: warnings for rows that don't have enough elements?
    // TODO: gender checking
    template<typename Range>
    bool Parse(const Range &text) {
        using namespace boost;
        using namespace std;
        // Construct the tokenizer
        typedef tokenizer<char_separator<char>,
                typename Range::const_iterator> tokenizer;
        char_separator<char> sep("\t", "\n", keep_empty_tokens);
        tokenizer tokens(text, sep);

        // Buffer to hold tokens
        // Use the 0-slot for unknown individuals
        vector<std::array<string, 6>> string_table(1);
        string_table.reserve(64);
        std::size_t k = 0;
        // Work through tokens and build vectors for each row
        for(auto tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
            if(*tok_iter == "\n") {
                k = 0;
                continue;
            } else if(k >= 6) {
                continue;
            } else if(k == 0) {
                string_table.emplace_back();
            }
            // Add token to the current row
            string_table.back()[k] = *tok_iter;
            k += 1;
        }
        // Build map of unique child names
        map<string, size_t> child_names;
        child_names.emplace("", 0);
        for(k = 1; k < string_table.size(); ++k) {
            bool success = child_names.emplace(string_table[k][1], k).second;
            // If child name is duplicate, erase it
            if(!success) {
                string_table[k][1].clear();
            }
        }
        // Parent->children relationships
        vector<vector<size_t>> children(string_table.size());
        for(k = 1; k < string_table.size(); ++k) {
            if(string_table[k][1].empty()) {
                continue;
            }
            // If parents are not known, use the 0-slot for them
            auto dad = child_names.find(string_table[k][2]);
            auto mom = child_names.find(string_table[k][3]);
            size_t ndad = (dad == child_names.end()) ? 0 : dad->second;
            size_t nmom = (mom == child_names.end()) ? 0 : mom->second;

            // If one parent is known and the other is unknown, fix it.
            // Use a \t to ensure no collision
            if(ndad == 0 && nmom != 0) {
                string par_name = "unknown_dad_of_" + string_table[k][1];
                string_table[k][2] = par_name;
                ndad = string_table.size();
                string_table.push_back({string_table[k][0], par_name,
                                        {}, {}, "male", {}
                                       });
                child_names.emplace(par_name, ndad);
                children.emplace_back();
            } else if(nmom == 0 && ndad != 0) {
                string par_name = "unknown_mom_of_" + string_table[k][1];
                string_table[k][3] = par_name;
                nmom = string_table.size();
                string_table.push_back({string_table[k][0], par_name,
                                        {}, {}, "female", {}
                                       });
                child_names.emplace(par_name, nmom);
                children.emplace_back();
            }

            // Push k onto dad and mom
            children[ndad].push_back(k);
            children[nmom].push_back(k);
        }

        // Add a dummy 0-th pedigree member to handle unknown individuals.
        // Use a breadth-first search starting at 0 to order pedigree
        names_.clear();
        names_.push_back("");
        vector<char> touched(string_table.size(), 0);
        touched[0] = 2;
        deque<size_t> visited{0};
        while(!visited.empty()) {
            size_t id = visited.front();
            visited.pop_front();
            for(auto a : children[id]) {
                // only add child to list after we have visited both parents
                if(touched[a] == 0) {
                    touched[a] = 1;
                } else if(touched[a] == 1) {
                    names_.push_back(string_table[a][1]);
                    visited.push_back(a);
                    touched[a] = 2;
                }
            }
        }
        // Construct table in the order of names_
        table_.clear();
        table_.emplace_back(Member{0, 0, 0, 0, Gender::Unknown, ""});
        for(auto && name : names_) {
            auto nrow = child_names[name];
            auto child = id(string_table[nrow][1]);


            if(child == 0) {
                continue;
            }
            table_.emplace_back(Member{
                0,
                child,
                id(string_table[nrow][2]),
                id(string_table[nrow][3]),
                ParseGender(string_table[nrow][4]),
                string_table[nrow][5]
            });
        }

        return true;
    }

protected:
    static Gender ParseGender(std::string &str) {
        static std::pair<std::string, Gender> keys[] = {
            {"0", Gender::Unknown},
            {"1", Gender::Male},
            {"2", Gender::Female},
            {"male", Gender::Male},
            {"female", Gender::Female},
            {"unknown", Gender::Unknown}
        };
        return dng::util::key_switch_tuple(str, keys, keys[0]).second;
    }

    NameContainer names_;
    MemberTable table_;
};

}
} // namespace dng::io

#endif // DNG_PEDIGREE_H
