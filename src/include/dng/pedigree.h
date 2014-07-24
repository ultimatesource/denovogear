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
#ifndef DNG_PEDIGREE_H
#define DNG_PEDIGREE_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

#include <boost/tokenizer.hpp>

namespace dng {

class Pedigree {
public:
	typedef boost::multi_index_container<std::string,
		boost::multi_index::indexed_by<
			boost::multi_index::random_access<>,
			boost::multi_index::ordered_unique<
			boost::multi_index::identity<std::string>>
		>> NameContainer;

	typedef std::vector<std::vector<std::string>> DataTable;
	
	Pedigree() {
	}
	
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
	
	// A reference to the data table
	const DataTable & table() const {
		return table_;
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
		
		// reset the pedigree string table
		table_.reserve(128);
		table_.resize(1);
		table_.back().clear();
		table_.back().reserve(6);
		
		// Work through tokens and build vectors for each row
		for (auto tok_iter = tokens.begin(); tok_iter != tokens.end();
			++tok_iter) {
			if(*tok_iter == "\n") {
				// Resize so that we have six elements
				table_.back().resize(6,"");
				// Push new row onto table
				table_.emplace_back();
				table_.back().reserve(6);
			} else {
				// Add token to the current row
				table_.back().push_back(*tok_iter);
			}
		}
		table_.back().resize(6,"");
		
		// Go through col 1 and pull out child names
		// Add a dummy 0-th pedigree member to handle
		// unknown individuals.
		names_.clear();
		names_.push_back("");
		for( auto &row : table_)
			names_.push_back(row[1]);
		
		return true;
	}
		
protected:
	NameContainer names_;
	DataTable table_;
};


} // namespace dng

#endif // DNG_PEDIGREE_H
