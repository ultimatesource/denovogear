/*
 * Copyright (c) 2015 Reed A. Cartwright
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

#include <dng/utility.h>

namespace dng {

class AlleleDepths {
public:
	typedef std::vector<int> data_t; 

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
	data_t::size_type num_libraries() const { return num_libraries_; }
	data_t::size_type num_nucleotides() const { return num_nucleotides_; }

protected:
	location_t location_;
	int8_t  type_;
	data_t data_;
	data_t::size_type num_nucleotides_;
	data_t::size_type num_libraries_;

};

}

#endif
