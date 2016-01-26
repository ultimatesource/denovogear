/*
 * Copyright (c) 2015-2016 Reed A. Cartwright
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
#ifndef DNG_DEPTHS_H
#define DNG_DEPTHS_H

#include <utility>

#include <dng/utility.h>

namespace dng {

typedef std::vector<int32_t> allele_depths_t; 

class AlleleDepths {
public:
	typedef allele_depths_t::size_type size_type;

	AlleleDepths(location_t location, int8_t type, size_type num_nucleotides, size_type num_libraries, allele_depths_t data) :
		location_(location), type_(type), num_nucleotides_(num_nucleotides), num_libraries_(num_libraries),
		data_(std::move(data))
	{
		assert(data_.size() == num_libraries_*num_nucleotides_);
	}

	AlleleDepths(location_t location, int8_t type, size_type num_nucleotides, size_type num_libraries) :
		location_(location), type_(type), num_nucleotides_(num_nucleotides), num_libraries_(num_libraries),
		data_(num_nucleotides*num_libraries, 0)
	{
	}

	int operator()(data_t::size_type nuc, data_t::size_type lib) const {
		// storage is nucleotide major
		assert(0 <= nuc && nuc < num_nucleotides_ && 0 <= lib && lib < num_libraries_);
		return data_[nuc*num_libraries_ + lib];
	}
	int& operator()(data_t::size_type nuc, data_t::size_type set) {
		// storage is nucleotide major
		assert(0 <= nuc && nuc < num_nucleotides_ && 0 <= lib && lib < num_libraries_);
		return data_[nuc*num_libraries_ + lib];
	}


	location_t location() const { return location_; }
	int8_t type() const { return type_; }
	size_type num_libraries() const { return num_libraries_; }
	size_type num_nucleotides() const { return num_nucleotides_; }

	void location(location_t location) { location_ = location; }
    void location(int target, int position) { location_ = make_location(target, position); }
	
	void type(int8_t type) { type_ = type; }

	void resize(size_type num_nucleotides, size_type num_libraries) {
		num_nucleotides_ = num_nucleotides;
		num_libraries_ = num_libraries;
		data_.resize(num_nucleotides*num_libraries);
	}
	void zero() {
		data_.assign(data_.size(), 0);
	}


protected:
	location_t location_;
	int8_t  type_;
	size_type num_nucleotides_;
	size_type num_libraries_;
	allele_depths_t data_;
};

} // namespace dng

#endif
