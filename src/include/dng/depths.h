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
#include <cstdint>
#include <unordered_map>

#include <boost/algorithm/string/case_conv.hpp>

#include <dng/utility.h>

namespace dng {
namespace pileup {

class AlleleDepths {
public:
    // information used to decode/encode AlleleDepths types
    struct type_info_t {
        char color; // binary id of the type
        char width; // the number of nucleotides in the type
        char label_upper[6]; // upper_case text version of the type
        char label_lower[6]; // lower_case text version of the type
        char reference; // the value of the reference base
        char indexes[4]; // the value of the alleles. Only use .width many. 
    };
    // information used to determine which genotypes are compatible with types
    struct type_info_gt_t {
        char color; // binary id of the type
        char width; // the number of nucleotides in the type
        char indexes[10]; // the value of the genotypes. Only use .width many.
    };

    static constexpr int type_info_table_length = 128;
    static const type_info_t type_info_table[128];
    static const type_info_gt_t type_info_gt_table[128];

    struct match_labels_t {
        std::unordered_map<std::string,int> tree;
        match_labels_t();
        int operator()(std::string str) const;
  
    };
    static match_labels_t MatchLabel;

    struct match_indexes_t {
        std::unordered_map<std::string,int> tree;
        match_indexes_t();

        int operator()(const std::string &rng) const;
    };
    static match_indexes_t MatchIndexes;

    static int8_t ColorDropN(int8_t color) { return color & 0x3F;}

    typedef std::vector<int32_t> data_t;
    typedef data_t::size_type size_type;

    AlleleDepths() : location_{0}, color_{0}, num_libraries_{0} { }

    AlleleDepths(location_t location, int8_t color, size_type num_lib, data_t data) :
        location_(location), color_(color), num_libraries_(num_lib),
        data_(std::move(data))
    {
        assert(data_.size() == num_libraries()*num_nucleotides());
    }

    AlleleDepths(location_t location, int8_t color, size_type num_lib) :
        location_(location), color_(color), num_libraries_(num_lib),
        data_(num_lib*type_info_table[color].width, 0)
    {
    }

    data_t::const_reference operator()(data_t::size_type nuc, data_t::size_type lib) const {
        // storage is nucleotide major
        assert(0 <= nuc && nuc < num_nucleotides() && 0 <= lib && lib < num_libraries());
        return data_[nuc*num_libraries_ + lib];
    }
    data_t::reference operator()(data_t::size_type nuc, data_t::size_type lib) {
        // storage is nucleotide major
        assert(0 <= nuc && nuc < num_nucleotides() && 0 <= lib && lib < num_libraries());
        return data_[nuc*num_libraries_ + lib];
    }

    location_t location() const { return location_; }
    void location(location_t location) { location_ = location; }
    void location(int contig, int position) { location_ = utility::make_location(contig, position); }

    int8_t color() const { return color_; }
    void color(int8_t color) { color_ = color; }
    const type_info_t& type_info() const { 
        assert(0 <= color_);
        return type_info_table[color_];
    }
    const type_info_gt_t& type_gt_info() const {
        assert(0 <= color_);
        return type_info_gt_table[color_];
    }

    const data_t& data() const { return data_; }
    data_t& data() { return data_; }
    size_type data_size() const { return data_.size(); };
    void data(data_t data) {
        assert(data.size() == data_.size());
        data_ = std::move(data);
    }
    void data_copy(const data_t& data) {
        assert(data.size() == data_.size());
        data_ = data;
    }
    void data_swap(data_t& data) {
        assert(data.size() == data_.size());
        data_.swap(data);
    }

    size_type num_libraries() const { return num_libraries_; }
    size_type num_nucleotides() const { return type_info().width; }
    std::pair<size_type,size_type> dimensions() const { return {num_nucleotides(), num_libraries()}; }
    void resize(int8_t color) {
        color_ = color;
        data_.resize(num_nucleotides()*num_libraries());
    }
    void resize(int8_t color, size_type num_lib) {
        color_ = color;
        num_libraries_ = num_lib;
        data_.resize(num_nucleotides()*num_libraries());
    }
    void zero() {
        data_.assign(data_.size(), 0);
    }
    void swap(AlleleDepths& other) {
        std::swap(location_, other.location_);
        std::swap(color_, other.color_);
        std::swap(num_libraries_, other.num_libraries_);
        data_.swap(other.data_);
    }

protected:
    location_t location_;
    size_type num_libraries_;
    data_t data_;
    int color_;
};

static_assert(sizeof(AlleleDepths::type_info_table) / sizeof(AlleleDepths::type_info_t) == AlleleDepths::type_info_table_length,
    "AlleleDepths::type_info_table does not have 128 elements.");
static_assert(sizeof(AlleleDepths::type_info_gt_table) / sizeof(AlleleDepths::type_info_gt_t) == AlleleDepths::type_info_table_length,
    "AlleleDepths::type_info_gt_table does not have 128 elements.");

inline
AlleleDepths::match_labels_t::match_labels_t() {
    for(int i=0; i<type_info_table_length; ++i) {
        auto & slot = AlleleDepths::type_info_table[i];
        tree.emplace(slot.label_upper, i);
    }
}

inline
int AlleleDepths::match_labels_t::operator()(std::string str) const {
    boost::to_upper(str);
    auto it = tree.find(str);
    return (it != tree.end()) ? it->second : -1;
}

inline
AlleleDepths::match_indexes_t::match_indexes_t() {
    // only add the first half of the table
    for(int i=0; i<type_info_table_length/2; ++i) {
        auto & slot = AlleleDepths::type_info_table[i];
        tree.emplace(std::string(&slot.indexes[0], &slot.indexes[slot.width]), i);
    }    
}

inline
int AlleleDepths::match_indexes_t::operator()(const std::string &rng) const {
    auto it = tree.find(rng);
    return (it != tree.end()) ? it->second : -1;
}

} // namespace pileup
} // namespace dng

#endif
