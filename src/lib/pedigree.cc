/*
 * Copyright (c) 2018 Reed A. Cartwright
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

#include <dng/pedigree.h>

namespace dng {

Pedigree Pedigree::parse_table(const std::vector<std::vector<std::string>> &table) {
   // Build Pedigree Object
    Pedigree ret;
    int row_num = 0;
    for(auto &&row : table) {
        Member member;
        row_num += 1;
        if(row.empty()) {
            continue;
        }
        if(row.size() < 5) {
            throw std::invalid_argument("Pedigree parsing failed. Row "
                + std::to_string(row_num) + " has "
                + std::to_string(row.size()) + " column(s) instead of 5 or more columns."
            );
        }
        // separate tags from member name
        {
            if(row[0][0] == '@') {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has empty child name."
                );                
            }
            auto tokens = utility::make_tokenizer_dropempty(row[0], "@", "");
            if(tokens.begin() == tokens.end()) {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has empty child name."
                );
            }
            auto it = tokens.begin();
            member.name = *(it++);
            member.tags.assign(it, tokens.end());
        }
        // process dad column
        {
            auto pos = row[1].find(':');
            member.dad = row[1].substr(0,pos);
            if(member.dad->empty()) {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has empty dad name."
                );
            }
            if(member.dad.get() == ".") {
                member.dad = boost::none;
            } else if(pos != std::string::npos) {
                char *str_end;
                member.dad_length = std::strtod(row[1].c_str()+pos+1, &str_end);
                if(pos+1 == row[1].size() || str_end != row[1].c_str()+row[1].size()) {
                    throw std::invalid_argument("Pedigree parsing failed. Row "
                        + std::to_string(row_num) + " has invalid dad length."
                    );
                }
            }
        }
        // process mom column
        {
            auto pos = row[2].find(':');
            member.mom = row[2].substr(0,pos);
            if(member.mom.get().empty()) {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has empty mom name."
                );
            }
            if(member.mom.get() == ".") {
                member.mom = boost::none;
            } else if(pos != std::string::npos) {
                char *str_end;
                member.mom_length = std::strtod(row[2].c_str()+pos+1, &str_end);
                if(pos+1 == row[2].size() || str_end != row[2].c_str()+row[2].size()) {
                    throw std::invalid_argument("Pedigree parsing failed. Row "
                        + std::to_string(row_num) + " has invalid mom length."
                    );
                }
            }
        }
        // process chromosmal sex column
        {
            member.sex = parse_sex(row[3]);
            if(member.sex == Sex::Unknown) {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has invalid sex."
                );
            }
        }
        // Process samples
        for(auto it = row.begin()+4; it != row.end();++it) {
             if(*it == "=") {
                member.samples.push_back(member.name);
            } else if(*it != ".") {
                member.samples.push_back(*it);
            }
        }
        ret.AddMember(member);
    }
    return ret;
}

} // namespace dng
