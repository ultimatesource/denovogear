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

#include <dng/utility.h>

namespace dng {
namespace pileup {

class AlleleDepths {
public:
    // information used to decode/encode AlleleDepths types
    struct type_info_t {
        char type;  // binary id of the type
        char width; // the number of nucleotides in the type
        char label_upper[6]; // upper_case text version of the type
        char label_lower[6]; // lower_case text version of the type
        char reference; // the value of the reference base
        char indexes[4]; // the value of the alleles.  Only use .width many.
    };
    static const type_info_t type_info_table[128];
    static int8_t LookupType(const std::vector<std::size_t> &indexes, bool ref_is_n);

    typedef std::vector<int32_t> data_t;
    typedef data_t::size_type size_type;

    AlleleDepths(location_t location, int8_t type, size_type num_lib, data_t data) :
        location_(location), type_(type), num_libraries_(num_lib),
        data_(std::move(data))
    {
        assert(data_.size() == num_libraries()*num_nucleotides());
    }

    AlleleDepths(location_t location, int8_t type, size_type num_lib) :
        location_(location), type_(type), num_libraries_(num_lib),
        data_(num_lib*type_info_table[type].width, 0)
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
    void location(int target, int position) { location_ = utility::make_location(target, position); }

    int8_t type() const { return type_; }
    void type(int8_t type) { type_ = type; }
    const type_info_t& type_info() const { return type_info_table[type_]; }


    const data_t& data() const { return data_; }
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
    void resize(size_type num_lib) {
        num_libraries_ = num_lib;
        data_.resize(num_nucleotides()*num_libraries());
    }
    void resize(int8_t type, size_type num_lib) {
        num_libraries_ = num_lib;
        type_ = type;
        data_.resize(num_nucleotides()*num_libraries());
    }
    void zero() {
        data_.assign(data_.size(), 0);
    }
    void swap(AlleleDepths& other) {
        std::swap(location_, other.location_);
        std::swap(type_, other.type_);
        std::swap(num_libraries_, other.num_libraries_);
        data_.swap(other.data_);
    }

protected:
    location_t location_;
    size_type num_libraries_;
    data_t data_;
    int8_t  type_;
};

} // namespace pileup
} // namespace dng

#endif
