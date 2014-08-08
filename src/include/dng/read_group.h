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
#include <boost/container/flat_set.hpp>

#include <unordered_set>

namespace dng {

class ReadGroups {
public:
	typedef boost::container::flat_set<std::string> StrSet;

	ReadGroups() { }
	
	template<typename InFiles>
	ReadGroups(InFiles &range) {
		Parse(range);
	}

	template<typename InFiles>
	void Parse(InFiles &range);

	template<typename Str>
	std::size_t group_index(const Str& id) const {
		auto it = groups_.find(id);
		return (it == groups_.end()) ? -1
			: static_cast<std::size_t>(it-groups_.begin());
	}
	
	inline const StrSet& groups() const { return groups_; }
	
	inline std::size_t library_index(std::size_t x) const {
		return library_ids_[x];
	}
	inline std::size_t sample_index(std::size_t x) const {
		return sample_ids_[x];
	}
	
	inline std::size_t num_groups() const { return groups_.size(); }
	inline std::size_t num_libraries() const { return libraries_.size(); }
	inline std::size_t num_samples() const { return samples_.size(); }

protected:

	StrSet groups_;
	StrSet libraries_;
	StrSet samples_;
	
	std::vector<std::size_t> library_ids_;
	std::vector<std::size_t> sample_ids_;
};

template<typename InFiles>
void ReadGroups::Parse(InFiles &range) {
	// structure to hold information on each read-group

	struct rg_t {
		std::string sample;
		std::string library;
	};
	std::map<std::string,rg_t> group_data;
	
	// iterate through each file in the range
	for(auto &f : range) {
		// get header text and continue on failre
		const char *text = f.header()->text;
		if(text == nullptr)
			continue;
		
		// enumerate over read groups
		for(text = strstr(text, "@RG\t"); text != nullptr;
		    text = strstr(text, "@RG\t"))
		{
			text += 4; // skip @RG\t
			// parse the @RG line as tab-separated key:value pairs
			// use a map to separate the parsing logic from the rg_t struct
			const char *k = text, *v=k,*p=k;
			std::map<std::string,std::string> tags;
			for(;*p != '\n' && *p != '\0';++p) {
				if(*p == ':')
					v = p+1; // make v point the the beginning of the value
					         // v-1 will point to the end of the key
				else if(*p == '\t') {
					if( k < v ) {
						// if we found a key:value pair, add it to the map
						// first one wins
						tags.emplace(std::string(k,v-1), std::string(v,p));
					}
					// move k and v to the next character
					k = v = p+1;
				}
			}
			// make sure the insert the last item if needed
			if( k < v )
				tags.emplace(std::string(k,v-1), std::string(v,p));
			
			// continue if this @RG does not have and ID
			auto it = tags.find("ID");
			if(it == tags.end() || it->second.empty())
				continue;
			// constuct the rg_t
			rg_t& val = group_data[it->second];
			// check to see if this is a duplicate
			if(!val.library.empty())
				continue;
			// library tag
			auto it2 = tags.find("LB");
			if(it2 != tags.end())
				val.library = std::move(it2->second);
			else
				val.library = it->second;
			// sample tag
			if((it2 = tags.find("SM")) != tags.end())
				val.sample = std::move(it2->second);
			else
				val.sample = it->second;
		}
	}
		
	// construct data vectors
	groups_.reserve(group_data.size());
	libraries_.reserve(group_data.size());
	samples_.reserve(group_data.size());	
	for(auto &a : group_data) {
		groups_.insert(a.first);
		libraries_.insert(a.second.library);
		samples_.insert(a.second.sample);
	}
	
	// connect groups to libraries and samples	
	library_ids_.resize(groups_.size(),-1);
	sample_ids_.resize(libraries_.size(),-1);
	size_t u=0;
	for(auto &a : group_data) {
		auto it = libraries_.find(a.second.library);
		assert(it != libraries_.end());
		library_ids_[u] = (it-libraries_.begin());
		it = samples_.find(a.second.sample);
		assert(it != samples_.end());
		sample_ids_[library_ids_[u]] = (it-samples_.begin());
		++u;
	}
}

};

#endif

