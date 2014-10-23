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
#ifndef DNG_IO_PED_H
#define DNG_IO_PED_H

#include <iostream>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

#include <boost/tokenizer.hpp>

namespace dng { namespace io {

class Pedigree {
public:
	typedef boost::multi_index_container<std::string,
		boost::multi_index::indexed_by<
			boost::multi_index::random_access<>,
			boost::multi_index::ordered_unique<
			boost::multi_index::identity<std::string>>
		>> NameContainer;

	typedef std::vector<std::vector<std::string>> DataTable;
	
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

	// Fetch the name of a member of the pedigree
	const std::string& name(std::size_t id) const {
		return names_[id];
	}
	// How many members, including the dummy, are in the pedigree.
	std::size_t member_count() const {
		return names_.size();
	}
	
	// Given the name of an individual return its id.
	// If the name is not valid, return the id of the 0-th
	// individual.
	std::size_t id(const std::string &name) const {
		auto it = names_.get<1>().find(name);
		if(it == names_.get<1>().end())
			return std::size_t(0);
		return names_.project<0>(it) - names_.begin();
	}

	std::vector<std::vector<std::string>> table() const {
		return {};
	}	
		
	// Parse a string-like object into a pedigree
	// TODO: maybe we can just keep pointers to the deliminators in memory
	// TODO: add comment support
	// TODO: warnings for rows that don't have enough elements?
	// TODO: gender checking
	// TODO: convert datatable to indexed relationships
	template<typename Range>
	bool Parse(const Range &text) {
		using namespace boost;
		using namespace std;
		// Construct the tokenizer
		typedef tokenizer<char_separator<char>,
			typename Range::const_iterator> tokenizer;
		char_separator<char> sep("\t", "\n", keep_empty_tokens);
		tokenizer tokens(text, sep);
		
		// Add a dummy 0-th pedigree member to handle
		// unknown individuals.
		names_.clear();
		row_ids_.clear();
		table_.clear()
		names_.push_back("");
		row_ids_.push_back(0);
		table_.push_back({0,0,0,0,Gender::Unknown,"");

		// Buffer to hold tokens
		vector<array<string,6>> string_table(1);
		string_table.reserve(64);
		std::size_t k = 0;
		// Work through tokens and build vectors for each row
		for (auto tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
			if(*tok_iter == "\n") {
				k = 0;
				continue;
			} else if(k >= 6) {
				continue;
			} else if(k == 0) {
				string_table.push_back();
			}
			// Add token to the current row
			string_table.back()[k] = *tok_iter;
			k += 1;
		}


		row_ids_.resize(names_.size(),0);

		// Process all parents and rewrite unknown names
		for(size_t u = 0; u < table_.size(); ++u) {
			auto & mrow = table_[u];
			// Point all unknown mom and dad entries to the unknown slot
			if(row_ids_[mrow.dad] == 0 && mrow.dad != 0)
				mrow.dad = 0;
			if(row_ids_[mrow.mom] == 0 && mrow.mom != 0)
				mrow.mom = 0;
			// If only one parent specified, create a dummy parent
			if(mrow.mom == 0 && mrow.dad != 0) {
				// construct a pseudo mom name using tab, which prevents collisions
				string mom_name = "dng:mom of\t" + names_[mrow.child];
				size_t id = names_.size();
				names_.push_back(mom_name);
				row_ids_.push_back(table_.size());
				table_.push_back({mrow.fam,id,0,0,Gender::Female,""});
				mrow.mom = id;
			}
			if(mrow.dad == 0 && mrow.mom != 0) {
				// construct a pseudo dad name using tab, which prevents collisions
				string dad_name = "dng:dad of\t" + names_[mrow.child];
				size_t id = names_.size();
				names_.push_back(dad_name);
				row_ids_.push_back(table_.size());
				table_.push_back({mrow.fam,id,0,0,Gender::Male,""});
				mrow.dad = id;
			}
		}
				

		return true;
	}
		
protected:
	inline bool AddRow(const std::vector<std::string>& row) {
		auto itc = names_.push_back(row[1]).first; // child
		size_t id = static_cast<size_t>(itc - names_.begin());
		// check to see if this child ID has been processed before
		if(id >= row_ids_.size()) {
			row_ids_.resize(id+1,0);
		} else if(row_ids_[id] != 0) {
			return false;
		}
		row_ids_[id] = table_.size();

		auto itf = names_.push_back(row[2]).first; // father
		auto itm = names_.push_back(row[3]).first; // mother
		table_.push_back({
			0,
			id,
			static_cast<size_t>(itf - names_.begin()),
			static_cast<size_t>(itm - names_.begin()),
			Gender::Unknown,
			row[5]
		});
		return true;
	}

	NameContainer names_;
	std::vector<std::size_t> row_ids_; // links child id to table slot
	MemberTable table_;
};

}} // namespace dng::io

#endif // DNG_PEDIGREE_H
