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
#ifndef DNG_PILEUP_H
#define DNG_PILEUP_H

#include <vector>
#include <limits>
#include <htslib/sam.h>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/metafunctions.hpp> 

#include <boost/intrusive/list.hpp>
#include <boost/noncopyable.hpp>

#include <dng/fileio.h>

namespace dng {

namespace detail {

using namespace boost::intrusive;

class SeqPool : boost::noncopyable {
public:
	struct node_t : public list_base_hook<> {
		bam1_t seq;  // sequence record;
		int32_t beg; // beginning of alignment
		int32_t end; // end of alignment
		int32_t pos; // current position of the pileup in this query/read
		int32_t rg;  // the read_group of the sequence	
			
		~node_t() {
			free(seq.data);
		}
	};
	typedef boost::intrusive::list<node_t> list_type;
	
	SeqPool(std::size_t sz=1024, std::size_t maxsz=std::numeric_limits<std::size_t>::max()) : 
		block_size_(sz/2), half_max_block_size_(maxsz/2)
	{
		if(block_size_ == 0)
			block_size_ = 1;
		
		store_.reserve(64);
		Expand();
	}
	
	node_t & Malloc() {
		// Check the inactive list first
		if(!inactive_.empty()) {
			node_t &n = inactive_.front();
			inactive_.pop_front();
			return n;
		}
		// Expand allocated space as needed
		if(next_ == end_)
			Expand();
		// Inplace allocation
		node_t *p = next_++;
		new (p) node_t();
		return *p;
	}
	
	void Free(node_t & n) {
		// Put the node on the list of free elements
		inactive_.push_front(n);
	}
	
	virtual ~SeqPool() {
		if(store_.empty())
			return;
		// enumerate over allocated blocks
		for(std::size_t i=0;i<store_.size()-1;++i) {
			for(std::size_t j=0;j<store_[i].second;++j) {
				// call destructor
				(store_[i].first+j)->~node_t();
			}
			// free memory
			std::return_temporary_buffer(store_[i].first);
		}
		// the last block is special
		for(node_t *p = store_.back().first;p!=next_;++p)
			p->~node_t();
		std::return_temporary_buffer(store_.back().first);
	}
	
protected:
	std::size_t block_size_;
	std::size_t half_max_block_size_;
	
	list_type inactive_;
	
private:
	void Expand() {
		// double block size until a limit is reached
		block_size_ = 2*std::min(block_size_, half_max_block_size_);
		// allocate space for more nodes
		auto a = std::get_temporary_buffer<node_t>(block_size_);
		if(a.first == nullptr)
			throw std::bad_alloc();
		next_ = a.first;
		end_  = next_ + a.second;
		store_.push_back(a);
	}
	
	std::vector<std::pair<node_t *,ptrdiff_t>> store_;
	node_t *next_, *end_;
};


} // namespace detail

class MPileup {
public:
	typedef std::vector<detail::SeqPool::list_type> data_type;
	typedef void (callback_type)(const data_type&,int,int);

	template<typename Input, typename Func>
	void operator()(Input &range, Func func);
	
protected:
	template<typename Input>
	int Advance(Input &range, data_type &data, int &target_id, int &ref_pos);

private:
	detail::SeqPool pool;

};


template<typename Input, typename Func>
void MPileup::operator()(Input &range, Func func) {
	using namespace std;
	using namespace fileio;
	
	// type erase our callback function
	function<callback_type> call_back(func);
	
	// data will hold our pileup information
	data_type data;
	int target_id = 0;
	int ref_pos = 0;
		
	// If there is a parsed region, use it.
	// TODO: check to see if the regions are all the same???
	int beg = 0, end = std::numeric_limits<int>::max();
	if(boost::begin(range)->iter() != nullptr) {
		beg = boost::begin(range)->iter()->beg;
		end = boost::begin(range)->iter()->end;
	}

	while(Advance(range,data,target_id, ref_pos) > 0) {
		if (pos < beg || pos >= end)
			continue; // out of range; skip
		//if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
		// Execute callback function
		call_back(data, target_id, pos);		
	}
}

template<typename Input>
int MPileup::Advance(Input &range, data_type &data, int &target_id, int &ref_pos) {
	/* Iterate through each element in range.
	     If the current position of range[i] equals the pileup minimum
	       read pileup queries from range[i] and advance position i
	       push reads onto a file-based buffer
	     After all streams have been advanced, update the minimum
	   Iterate through the buffer of each file
	     push reads onto the approrpraite read-group in data
	     update positions
	*/
}



typedef void (mpileup_func_t)(int,int,const std::vector<int>&,
	const std::vector<const bam_pileup1_t *>&);

template<typename Input, typename Func>
void mpileup(Input &range, bam_plp_auto_f auto_func, Func func) {
	using namespace std;
	using namespace fileio;
		
	// type erase our callback function
	function<mpileup_func_t> call_back(func);
	
	// bam_mplp_init requires an array of pointers
	std::size_t num_files = boost::size(range);
	vector<typename boost::range_pointer<Input>::type>
		pointers(num_files, nullptr);
	{ std::size_t i=0;
	for(auto a=boost::begin(range); a != boost::end(range); ++a) {
		pointers[i] = &(*a) + i;
		i++;
	}}
	
	// If there is a parsed region, use it.
	// TODO: check to see if the regions are all the same???
	int beg = 0, end = std::numeric_limits<int>::max();
	if(boost::begin(range)->iter() != nullptr) {
		beg = boost::begin(range)->iter()->beg;
		end = boost::begin(range)->iter()->end;
	}
	
	// initialize the pileup iterator
	bam_mplp_t mpileup_iter = bam_mplp_init(num_files, auto_func,
		reinterpret_cast<void**>(pointers.data()));
	
	vector<int> counts(num_files, 0);
	vector<const bam_pileup1_t *> records(num_files, nullptr);
	int target_id, pos;
	
	while(bam_mplp_auto(mpileup_iter, &target_id, &pos, counts.data(),records.data()) > 0) {
		if (pos < beg || pos >= end) continue; // out of range; skip
		//if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
		call_back(target_id, pos, counts, records);
	}
	bam_mplp_destroy(mpileup_iter);
}

} //namespace dng

#endif //DNG_PILEUP_H

