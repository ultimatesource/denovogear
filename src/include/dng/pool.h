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
#ifndef DNG_POOL_H
#define DNG_POOL_H

#include <boost/intrusive/list.hpp>
#include <boost/noncopyable.hpp>

namespace dng {

namespace detail {
#ifdef NDEBUG
typedef boost::intrusive::link_mode<boost::intrusive::normal_link>
node_link_mode_t;
#else
typedef boost::intrusive::link_mode<boost::intrusive::safe_link>
node_link_mode_t;
#endif
}

typedef boost::intrusive::list_base_hook<detail::node_link_mode_t> PoolNode;

template<typename Node>
class IntrusivePool : boost::noncopyable {
public:
    typedef Node node_type;
    typedef boost::intrusive::list<node_type> list_type;

    IntrusivePool(std::size_t sz = 1024, std::size_t maxsz
                  = std::numeric_limits<std::size_t>::max()) :
        block_size_(sz / 2), half_max_block_size_(maxsz / 2),
        count_(0) {
        if(block_size_ == 0) {
            block_size_ = 1;
        }

        store_.reserve(64);
        Expand();
    }

    node_type *Malloc() {
        // Check the inactive list first
        if(!inactive_.empty()) {
            node_type *p = &inactive_.front();
            inactive_.pop_front();
            return p;
        }
        // Expand allocated space as needed
        if(next_ == end_) {
            Expand();
        }
        // Inplace allocation
        node_type *p = next_++;
        new(p) node_type();
        return p;
    }

    void Free(node_type *p) noexcept {
        // Put the node on the list of free elements
        inactive_.push_front(*p);
    }

    virtual ~IntrusivePool() {
        if(store_.empty()) {
            return;
        }
        // clear inactive as needed
        if(list_type::safemode_or_autounlink) {
            inactive_.clear();
        }
        // enumerate over allocated blocks
        for(std::size_t i = 0; i < store_.size() - 1; ++i) {
            for(std::size_t j = 0; j < store_[i].second; ++j) {
                // call destructor
                (store_[i].first + j)->~node_type();
            }
            // free memory
            std::return_temporary_buffer(store_[i].first);
        }
        // the last block is special
        for(node_type *p = store_.back().first; p != next_; ++p) {
            p->~node_type();
        }
        std::return_temporary_buffer(store_.back().first);
    }

protected:
    std::size_t block_size_;
    std::size_t half_max_block_size_;

    std::size_t count_;

    list_type inactive_;

private:
    void Expand() {
        // double block size until a limit is reached
        block_size_ = 2 * std::min(block_size_, half_max_block_size_);
        // allocate space for more nodes
        auto a = std::get_temporary_buffer<node_type>(block_size_);
        if(a.first == nullptr) {
            throw std::bad_alloc();
        }
        // set block size to actual size allocated
        block_size_ = a.second;
        // storage information
        next_ = a.first;
        end_  = next_ + a.second;
        store_.push_back(a);
    }

    std::vector<std::pair<node_type *, ptrdiff_t>> store_;
    node_type *next_, * end_;
};

} // namespace dng

#endif //DNG_POOL_H
