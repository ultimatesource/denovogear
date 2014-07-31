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
#include <unordered_map>
#include <htslib/sam.h>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/metafunctions.hpp> 

#include <boost/intrusive/list.hpp>
#include <boost/noncopyable.hpp>

#include <dng/fileio.h>
#include <dng/cigar.h>

namespace dng {

inline uint64_t make_location(int t, int p) {
	return (static_cast<uint64_t>(t) << 32) | p;
}
inline int location_to_target(uint64_t u) {
	return static_cast<int>(u >> 32);
}
inline int location_to_position(uint64_t u) {
	return static_cast<int>(u & 0x7FFFFFFF);
}

namespace detail {

/*
 @abstract Structure for core alignment information.
 @field  tid     chromosome ID, defined by bam_hdr_t
 @field  pos     0-based leftmost coordinate
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_qname length of the query name
 @field  flag    bitwise flag
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template

typedef struct {
	int32_t tid;
	int32_t pos;
	uint32_t bin:16, qual:8, l_qname:8;
	uint32_t flag:16, n_cigar:16;
	int32_t l_qseq;
	int32_t mtid;
	int32_t mpos;
	int32_t isize;
} bam1_core_t;

 @abstract Structure for one alignment.
 @field  core       core information about the alignment
 @field  l_data     current length of bam1_t::data
 @field  m_data     maximum length of bam1_t::data
 @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
 
 @discussion Notes:
 
 1. qname is zero tailing and core.l_qname includes the tailing '\0'.
 2. l_qseq is calculated from the total length of an alignment block
 on reading or from CIGAR.
 3. cigar data is encoded 4 bytes per CIGAR operation.
 4. seq is nybble-encoded according to bam_nt16_table.

typedef struct {
	bam1_core_t core;
	int l_data, m_data;
	uint8_t *data;
#ifndef BAM_NO_ID
	uint64_t id;
#endif
} bam1_t;
*/

using namespace boost::intrusive;

class SeqPool : boost::noncopyable {
public:
	struct node_t : public list_base_hook<> {
		bam1_t seq;   // sequence record : [beg,end)
		uint64_t beg; // 0-based left-most edge 
		uint64_t end; // 0-based right-most edge
		uint64_t pos; // current position of the pileup in this query/read
		bool is_missing; // is there no base call at the pileup position?
				
		int32_t target() const {
			return seq.core.tid;
		}
		cigar_t cigar() const {
			return std::make_pair(seq.core.n_cigar,bam_get_cigar(&seq));
		}
		
		~node_t() {
			free(seq.data);
		}
	};
	typedef boost::intrusive::list<node_t> list_type;
	typedef list_type::value_type node_type;
	
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
		inactive_.clear(); // TODO: determine best way to avoid this?
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
		// set block size to actual size allocated
		block_size_ = a.second;
		// storage information
		next_ = a.first;
		end_  = next_ + a.second;
		store_.push_back(a);
	}
	
	std::vector<std::pair<node_t *,ptrdiff_t>> store_;
	node_t *next_, *end_;
};

template<typename InFile>
class BamScan {
public:
	typedef SeqPool::list_type list_type;
	typedef SeqPool::node_type node_type;
	
	explicit BamScan(InFile& in) : in_(in), next_loc_(0) {
		
	}
	
	list_type operator()(uint64_t target_loc, SeqPool &pool) {
		// Reads from in_ until it encounters the first read
		// that is right of pos.
		// TODO: check for proper read ordering???
		while(next_loc_ <= target_loc) {
			node_type &n = pool.Malloc();
			do {
				// Try to grab a read, if not return current buffer
				if(in_(n.seq) < 0) {
					pool.Free(n);
					next_loc_ = std::numeric_limits<uint64_t>::max();
					return std::move(buffer_);
				}
				n.beg = make_location(n.seq.core.tid, n.seq.core.pos);
				// update right-most position in the read
				n.end = n.beg + cigar::target_length(n.cigar());
			} while( n.end <= target_loc);
			
			// update the left-most position of the most recently read read.
			next_loc_ = n.beg;
			// save read
			buffer_.push_back(n);	
		}
		// Return all but the last read.
		if(buffer_.size() <= 1)
			return list_type();
		list_type ret;
		ret.splice(ret.end(),buffer_,buffer_.begin(),--buffer_.end());
		return ret;
	}
	
	uint64_t next_loc() const { return next_loc_; }
	
private:
	uint64_t next_loc_;	
	InFile & in_;
	list_type buffer_;
};

} // namespace detail

class MPileup {
public:
	typedef detail::SeqPool::list_type list_type;
	typedef detail::SeqPool::node_type node_type;
	typedef std::vector<list_type> data_type;
	typedef void (callback_type)(const data_type&, uint64_t);

	template<typename InFiles, typename Func>
	void operator()(InFiles &range, Func func);
	
protected:
	template<typename Scanners>
	bool Advance(Scanners &range, data_type &data, uint64_t &target_loc,
		uint64_t &fast_forward_loc);

private:
	detail::SeqPool pool_;
	
	std::unordered_map<std::string,int> read_groups_;
};

template<typename InFiles, typename Func>
void MPileup::operator()(InFiles &range, Func func) {
	using namespace std;
	using namespace fileio;
	
	// encapsulate input files into scanners
	vector<detail::BamScan<typename boost::range_value<InFiles>::type>> 
		scanners(boost::begin(range),boost::end(range));
	
	// type erase our callback function
	function<callback_type> call_back(func);
	
	// data will hold our pileup information
	data_type data;
	uint64_t current_loc = 0;
	uint64_t fast_forward_loc = 0;
		
	// If there is a parsed region, use it.
	// TODO: check to see if the regions are all the same???
	uint64_t beg_loc = 0, end_loc = std::numeric_limits<uint64_t>::max();
	if(boost::begin(range)->iter() != nullptr) {
		beg_loc = make_location(boost::begin(range)->iter()->tid,
			boost::begin(range)->iter()->beg);
		end_loc = make_location(boost::begin(range)->iter()->tid,
			boost::begin(range)->iter()->end);
	}
	
	while(Advance(scanners,data,current_loc,fast_forward_loc)) {
		if (current_loc < beg_loc || current_loc >= end_loc)
			continue; // out of range; skip
		//if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
		// Execute callback function
		call_back(data, current_loc);
	}
}

// Advance starts a pileup procedure at pos.
// If there are no reads at pos, it forwards to the next pos and updates
// current_pos

template<typename Scanners>
bool MPileup::Advance(Scanners &range, data_type &data, uint64_t &target_loc,
	uint64_t &fast_forward_loc) {
	using namespace std;
	// enumerate through existing data set, updating location of reads and
	// purge reads that have expired
	bool fast_forward = true;
	for(auto &d : data) {
		auto it = d.begin();
		while(it != d.end()) {
			if(it->end <= target_loc) {
				d.erase(it++);
				continue;
			}
			uint64_t q = cigar::target_to_query(target_loc,it->beg,it->cigar());
			it->pos = cigar::query_pos(q);
			it->is_missing = cigar::query_del(q);
			++it;
		}
		// we want to fast_forward if there is nothing to output
		fast_forward = (fast_forward && d.empty());
	}
	if(fast_forward)
		target_loc = fast_forward_loc;
		
	if(target_loc >= fast_forward_loc) {
		uint64_t next_loc = numeric_limits<uint64_t>::max();
		for(auto &scanner : range) {
			// Scan reads from file until target_loc
			list_type new_reads = scanner(target_loc,pool_);
			// Update the minimum position of the next read
			next_loc = std::min(next_loc,scanner.next_loc());
			// process read_groups
			while(!new_reads.empty()) {
				node_type &n = new_reads.front();
				new_reads.pop_front();
				uint8_t *rg = bam_aux_get(&n.seq, "RG");
				auto it = read_groups_.find(reinterpret_cast<const char*>(rg+1));
				if(it == read_groups_.end()) {
					pool_.Free(n); // drop unknown RG's
					continue;
				}
				// process cigar string
				uint64_t q = cigar::target_to_query(target_loc,n.beg,n.cigar());
				n.pos = cigar::query_pos(q);
				n.is_missing = cigar::query_del(q);

				// push read onto correct RG
				data[it->second].push_back(n);
				fast_forward = false;
			}
		}
		fast_forward_loc = next_loc;
	}
	return (!fast_forward);
}

/*
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
*/

} //namespace dng

#endif //DNG_PILEUP_H

