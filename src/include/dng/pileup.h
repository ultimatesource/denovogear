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
		bam1_t seq;  // sequence record;
		uint64_t beg;  // 0-based left-most coordinate
		uint64_t end; // 0-based right-most edge
		uint64_t pos; // current position of the pileup in this query/read
				
		int32_t target() const {
			return n.seq.core.tid;
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

inline uint64_t make_position(int t, int p) {
	return (static_cast<uint64_t>(t) << 32) | p;
}

template<typename InFile>
class BamScan {
public:
	typedef SeqPool::list_type list_type;
	typedef SeqPool::node_type node_type;
	
	explicit BamScan(InFile& in) : in_(in), pos_(0) {
		
	}
	
	list_type&& operator()(uint64_t target_pos, SeqPool &pool) {
		// Reads from in_ until it encounters the first read
		// that is right of pos.
		// TODO: check for proper read ordering???
		while(next_pos_ <= target_pos) {
			node_type &n = pool.Malloc();
			for(;;) {
				// Try to grab a read, if not return current buffer
				if(in_(n.seq) < 0) {
					pool.Free(n);
					next_pos_ = std::numeric_limits<uint64_t>::max();
					return std::move(buffer_);
				}
				n.beg = make_position(n.seq.core.tid, n.seq.core.pos)
				// update right-most position in the read
				n.right_pos = n.beg + cigar::target_length(n.seq.core.n_cigar,
					bam_get_cigar(&n.seq));
				if(n.right_pos > target_pos)
					break;
			}
			// update the left-most position of the most recently read read.
			next_pos_ = n.beg;
			// save read
			buffer_.push_back(n);	
		}
		// Return all but the last read.
		if(buffer_.size() <= 1)
			return std::move(list_type());
		list_type ret;
		ret.splice(ret.end(),buffer_,buffer_.begin(),--buffer_.end());
		return std::move(ret);
	}
	
	uint64_t next_pos() const { return next_pos_; }
	
private:
	uint64_t next_pos_;	
	InFile & in_;
	list_type buffer_;
};

void process_cigar(SeqPool::node_type &n);

} // namespace detail

class MPileup {
public:
	typedef std::vector<detail::SeqPool::list_type> data_type;
	typedef void (callback_type)(const data_type&,int,int);

	template<typename InFiles, typename Func>
	void operator()(InFiles &range, Func func);
	
protected:
	template<typename Scanners>
	int Advance(Scanners &range, data_type &data, uint64_t &target_pos),
		uint64_t &fast_forward_pos);

private:
	detail::SeqPool pool_;
	
	std::unordered_map<std::string,int> read_groups_;
};

template<typename InFiles, typename Func>
void MPileup::operator()(InFiles &range, Func func) {
	using namespace std;
	using namespace fileio;
	
	// encapsulate input files into scanners
	vector<BamScan> scanners(boost::begin(range),boost::end(range));
	
	// type erase our callback function
	function<callback_type> call_back(func);
	
	// data will hold our pileup information
	data_type data;
	uint32_t current_pos = 0;
		
	// If there is a parsed region, use it.
	// TODO: check to see if the regions are all the same???
	int beg = 0, end = std::numeric_limits<int>::max();
	if(boost::begin(range)->iter() != nullptr) {
		beg = boost::begin(range)->iter()->beg;
		end = boost::begin(range)->iter()->end;
	}
	
	while(Advance(scanners,data,current_pos) > 0) {
		if (pos < beg || pos >= end)
			continue; // out of range; skip
		//if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
		// Execute callback function
		call_back(data, target_id, pos);		
	}
}

// Advance starts a pileup procedure at pos.
// If there are no reads at pos, it forwards to the next pos and updates
// current_pos

template<typename Scanners>
int MPileup::Advance(Scanners &range, data_type &data, uint64_t &target_pos
	uint64_t &fast_forward_pos) {
	/* Iterate through each element in range.
	     If the current position of range[i] equals the pileup minimum
	       read pileup queries from range[i] and advance position i
	       push reads onto a file-based buffer
	     After all streams have been advanced, update the minimum
	   Iterate through the buffer of each file
	     push reads onto the approrpraite read-group in data
	     update positions
	*/
	using namespace std;
	// enumerate through existing data set, updating position of reads
	// purge reads that have expired
	bool fast_forward = true;
	for(auto d : data) {
		auto it = d.begin()
		while(it != d.end()) {
			if(it->end <= target_pos) {
				d.erase(it++);
				continue;
			}
			++it;
		}
		// we want to fast_forward if there is nothing to output
		fast_forward = fast_forward && d.empty();
	}
	if(fast_forward)
		target_pos = fast_forward_pos;
		
	if(target_pos >= fast_forward_pos) {
		uint64_t next_pos = numeric_limits<uint64_t>::max();
		for(auto scanner : range) {
			// Scan reads from file until target_pos
			list_type new_reads = scanner(target_pos,pool_);
			next_pos = std::min(next_pos,scanner.next_pos());
			// process read_groups
			while(!new_reads.empty()) {
				node_type &n = new_reads.front();
				new_reads.pop_front();
				uint8_t *rg = bam_aux_get(&n.seq, "RG");
				auto it = read_groups_.find(reinterpret_cast<const char*>(rg+1));
				if(it == read_groups_.end())
					pool_.Free(n); // drop unknown RG's
				else
					data[it.second].push_back(n); // push read onto correct RG
				// process cigar string
				
				
			}
		}
		fast_forward_pos = next_pos;
	}
	
	
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

