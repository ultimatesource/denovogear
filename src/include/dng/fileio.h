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

#include <htslib/hts.h>
#include <htslib/sam.h>

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

class SamFile : boost::noncopyable {
public:
	SamFile(const char *file, const char *region, int min_mapQ = 0, int min_len = 0) :
			fp_(nullptr), hdr_(nullptr), iter_(nullptr),
	    	min_mapQ_(min_mapQ), min_len_(min_len) {
	    fp_ = sam_open(file, "r");
	    if(fp_ == nullptr)
	    	throw std::runtime_error("unable to open file '" + std::string(file) + "'.");
	    hdr_ = sam_hdr_read(fp_);
	    if(region != nullptr && region[0] != '\0') {
	    	hts_idx_t *idx = sam_index_load(fp_, file);
	    	if(idx == nullptr)
	    		throw std::runtime_error("unable to load index for '" + std::string(file) + "'.");
	    	iter_ = sam_itr_querys(idx, hdr_, region);
	    	hts_idx_destroy(idx);
	    	if(iter_ == nullptr) {
	    		throw std::runtime_error("unable to parse region '" + std::string(region) + "'.");
	    	}
	    }
	}
	SamFile(SamFile &&other) : fp_(other.fp_), hdr_(other.hdr_),
		iter_(other.iter_)
	{
		other.fp_ = nullptr;
		other.hdr_ = nullptr;
		other.iter_ = nullptr;		
	}
	
	// move-and-swap
	SamFile& operator=(SamFile &other) {
		std::swap(fp_,other.fp_);
		std::swap(hdr_,other.hdr_);
		std::swap(iter_,other.iter_);
		return *this;
	}
	
	int operator()(bam1_t &b) {
		int ret;
		for(;;) {
		    ret = (iter_) ? sam_itr_next(fp_, iter_, &b)
		                  : sam_read1(fp_, hdr_, &b);
		    if (ret < 0)
		    	break;
		    if (b.core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))
		    	continue;
		    if (b.core.qual < min_mapQ_)
		    	continue;
		    if (min_len_ && bam_cigar2qlen(b.core.n_cigar, bam_get_cigar(&b)) < min_len_)
		    	continue;
		    break;
		}
		return ret;
	}

	// cleanup as needed
	virtual ~SamFile() {
		if(iter_ != nullptr)
			hts_itr_destroy(iter_);
		if(hdr_ != nullptr)
			bam_hdr_destroy(hdr_);
		if(fp_ != nullptr)
			sam_close(fp_);
	}
	
	const bam_hdr_t * header() const { return hdr_; }
	const hts_itr_t * iter() const { return iter_; }
	
protected:
	samFile *fp_;     // the file handle
	bam_hdr_t *hdr_;  // the file header
	hts_itr_t *iter_; // NULL if a region not specified
	int min_mapQ_, min_len_; // mapQ filter; length filter
};

} // namespace fileio

using fileio::ParsedList;
using fileio::SamFile;

} // namespace dng

#endif // DNG_FILEIO_H

