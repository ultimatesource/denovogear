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

	struct rg_t {
		std::string id;
		std::string library;
		std::string sample;
	};

	typedef boost::multi_index_container<rg_t, indexed_by<
		ordered_unique<tag<rg::id>, member<rg_t,std::string,&rg_t::id>>,
		ordered_non_unique<tag<rg::lb>, member<rg_t,std::string,&rg_t::library>>,
		ordered_non_unique<tag<rg::sm>, member<rg_t,std::string,&rg_t::sample>>
		>> DataBase;
}

 
class ReadGroups {
public:
	typedef boost::container::flat_set<std::string> StrSet;
	typedef detail::rg_t ReadGroup;
	typedef detail::DataBase DataBase;

	//ReadGroups() { }
	
	///template<typename InFiles>
	//ReadGroups(InFiles &range) {
	//	Parse(range);
	//}

	// Parse a single vcf file into read groups
	//void ParseVCF(const char *fname);


	
	// Parse a list of BAM/SAM files into readgroups
	template<typename InFiles>
	void Parse(InFiles &range);

	template<typename InFile>
	void Parse(InFile *fname);
	
	/*
	template<>
	  void Parse<const char>(const char *fname)
	  {

	  }
	*/
	
	inline const StrSet& groups() const { return groups_; }
	inline const StrSet& libraries() const { return libraries_; }
	inline const StrSet& samples() const { return samples_; }
	inline const DataBase& data() const { return data_; }

	inline std::size_t library_from_id(std::size_t n) const {
		return library_from_id_[n];
	}


protected:
	DataBase data_;

	StrSet groups_;
	StrSet libraries_;
	StrSet samples_;
	
	std::vector<std::size_t> library_from_id_;
};

namespace rg {
inline std::size_t index(const ReadGroups::StrSet& set, const std::string& query) {
	auto it = set.find(query);
	return (it != set.end()) ? static_cast<std::size_t>(it-set.begin()) : -1;
}
}


template<typename InFiles>
void ReadGroups::Parse(InFiles &range) {
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
			
			// continue to next @RG if this @RG does not have an ID
			auto it = tags.find("ID");
			if(it == tags.end() || it->second.empty())
				continue;
			// check to see if this is a duplicate
			// TODO: Raise warning/error if this is a collision of different sm/lb data
			if(data_.find(it->second) != data_.end())
				continue;
			ReadGroup val{it->second};

			// sample tag
			if((it = tags.find("SM")) != tags.end())
				val.sample = std::move(it->second);
			else
				val.sample = val.id;
			// library tag
			it = tags.find("LB");
			if(it != tags.end()) {
				val.library = std::move(it->second);
				// Prevent collision of LB tags from different samples by appending sample tag
				val.library += "\t" + val.sample;
			} else
				val.library = val.id + "\t" + val.sample;
			data_.insert(std::move(val));
		}
	}
	
	// construct data vectors
	groups_.reserve(data_.size());
	libraries_.reserve(data_.size());
	samples_.reserve(data_.size());	
	for(auto &a : data_) {
		groups_.insert(a.id);
		libraries_.insert(a.library);
		samples_.insert(a.sample);
	}
	
	// connect groups to libraries and samples	
	library_from_id_.resize(groups_.size());
	for(auto &a : data_) {
		library_from_id_[rg::index(groups_,a.id)] = rg::index(libraries_,a.library);

	}
}

template<typename InFFile>
void ReadGroups::Parse(InFFile *fname)
{

}
 
}

#endif

