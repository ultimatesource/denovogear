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
#ifndef DNG_FILEIO_H
#define DNG_FILEIO_H

#include <cstdlib>
#include <stdexcept>

#include <dng/hts/bam.h>

#include <boost/noncopyable.hpp>

namespace dng {
namespace fileio {
// NOTE: Most of fileio is wraps samtools and htslib.  Some code has been copied
//       from those libraries to support the linkage.

class ParsedList {
public:
	static const bool kFile = true;
	static const bool kString = false;
	
	ParsedList(const char *str, bool is_file=false) {
		list_ = hts_readlist(str, is_file ? 1 : 0, &length_);
	}
	
	bool Empty() const {
		return list_ == nullptr;
	}

	char * operator[](std::size_t k) {
		assert(k < length_); // check to see if k is valid
		return list_[k];
	}
	
	std::size_t Size() const {
		return length_;
	}	
	
	const char * operator[](std::size_t k) const {
		assert(k < length_); // check to see if k is valid
		return list_[k];
	}

	virtual ~ParsedList() {
		for(int i=0;i<length_;++i)
			free(list_[i]);
		if(length_ > 0)
			free(list_);
	}
private:
	char **list_;
	int length_;

};

} // namespace fileio

using fileio::ParsedList;

} // namespace dng

#endif // DNG_FILEIO_H

