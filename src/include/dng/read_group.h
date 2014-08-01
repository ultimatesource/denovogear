/*
 * Copyright (c) 2014 Reed A. Cartwright
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

#include <unordered_set>

namespace dng {

class ReadGroups {
public:
	template<typename InFiles>
	void Parse(InFiles &range);

	template<typename str>
	std::size_t Index(const str& id) {
		auto it = boost::lower_bound(ids_,id);
		return (it == ids_.end()) ? -1
			: static_cast<std::size_t>(it-boost::begin(ids_));
	}
	
	const std::vector<std::string>& all() const { return ids_; }

protected:
	std::vector<std::string> ids_;
};

template<typename InFiles>
void ReadGroups::Parse(InFiles &range) {
	ids_.clear();
	for(auto &f : range) {
		const char *text = f.header()->text;
		if(text == nullptr)
			continue;
		text = strstr(text, "@RG\t");
		while(text != nullptr) {
			text += 4;
			const char *k = text, *v=k,*p=k;
			std::map<std::string,std::string> tags;			
			for(;*p != '\n' && *p != '\0';++p) {
				if(*p == ':')
					v = p+1;
				else if(*p == '\t') {
					if( k < v ) {
						// first one wins
						tags.emplace(std::string(k,v-1), std::string(v,p));
					}
					k = v = p+1;
				}
			}
			if( k < v )
				tags.emplace(std::string(k,v-1), std::string(v,p));
			auto it_id = tags.find("ID");
			if(it_id != tags.end()) {
				ids_.push_back(it_id->second);
			}
			
			// find next readgroup
			text = strstr(text, "@RG\t");
		}
	}
	boost::erase(ids_, boost::unique<boost::return_found_end>(boost::sort(ids_)));
}

};

#endif

