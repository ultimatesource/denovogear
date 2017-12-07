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
#include <boost/functional/hash.hpp>

#include <dng/utility.h>

namespace dng {
namespace pileup {

struct depth_t {
    int32_t counts[4]{0,0,0,0};
};

using RawDepths = std::vector<depth_t>;

namespace detail {
template <typename Container>
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};    

};

class AlleleDepths {
public:
    // information used to decode/encode AlleleDepths types
    struct type_info_t {
        char color; // binary id of the type
        char width; // the number of nucleotides in the type
        char label_upper[6]; // upper_case text version of the type
        char label_lower[6]; // lower_case text version of the type
        char label_htslib[10]; // String used for htslib's bcf_update_alleles_str
        char reference; // the value of the reference base
        char indexes[4]; // the value of the alleles. Only use .width many. 
    };
    // information used to determine which genotypes are compatible with types
    struct type_info_gt_t {
        char color; // binary id of the type
        char width; // the number of nucleotides in the type
        char indexes[10]; // the value of the genotypes. Only use .width many.
        char haploid_indexes[4]; // location of homozygotes that are compatible with type
    };

    static constexpr int type_info_table_length = 128;
    static const type_info_t type_info_table[128];
    static const type_info_gt_t type_info_gt_table[128];

    static const char hash_to_color[256];

    static const int alleles_diploid[10][2];
    static const int encoded_alleles_diploid_unphased[10][2];
    static const int encoded_alleles_haploid[4][2];

    struct match_labels_t {
        std::unordered_map<std::string,int> tree;
        match_labels_t();
        int operator()(std::string str) const;
  
    };
    static match_labels_t MatchLabel;

    struct match_indexes_t {
        using key_t = std::vector<char>;
        std::unordered_map<key_t,int, detail::container_hash<key_t>> tree;
        match_indexes_t();

        int operator()(const key_t &rng) const;
    };
    static match_indexes_t MatchIndexes;

    static int MatchAlleles(char *str);

    static int8_t ColorDropN(int8_t color) { return color & 0x3F;}

    using data_t = std::vector<int32_t>;
    using size_type = data_t::size_type;

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

    data_t::reference operator()(data_t::size_type lib, data_t::size_type nuc) {
        // storage is library major
        assert(data_.size() == num_nucleotides()*num_libraries());
        assert(0 <= nuc && nuc < num_nucleotides() && 0 <= lib && lib < num_libraries());
        return data_[lib*num_nucleotides()+nuc];
    }

    data_t::const_reference operator()(data_t::size_type lib, data_t::size_type nuc) const {
        // storage is library major
        assert(data_.size() == num_nucleotides()*num_libraries());
        assert(0 <= nuc && nuc < num_nucleotides() && 0 <= lib && lib < num_libraries());
        return data_[lib*num_nucleotides()+nuc];
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
    data_t& data() { 
        assert(data_.size() == num_nucleotides()*num_libraries());
        return data_;
    }
    size_type data_size() const { return data_.size(); };

    size_type num_libraries() const { return num_libraries_; }
    size_type num_nucleotides() const { return type_info().width; }
    std::pair<size_type,size_type> dimensions() const { return {num_nucleotides(), num_libraries()}; }
    
    void Resize(int8_t color) {
        color_ = color;
        data_.resize(num_nucleotides()*num_libraries());
    }
    void Resize(int8_t color, size_type num_lib) {
        color_ = color;
        num_libraries_ = num_lib;
        data_.resize(num_nucleotides()*num_libraries());
    }
    void Zero() {
        assert(data_.size() == num_nucleotides()*num_libraries());
        data_.assign(data_.size(), 0);
    }
    size_type TotalDepth() const {
        assert(data_.size() == num_nucleotides()*num_libraries());
        size_type dp = 0;
        for(auto &&a : data_) {
            dp += a;
        }
        return dp;
    }

    void swap(AlleleDepths& other) {
        std::swap(location_, other.location_);
        std::swap(color_, other.color_);
        std::swap(num_libraries_, other.num_libraries_);
        data_.swap(other.data_);
    }

protected:
    location_t location_;
    int color_;
    size_type num_libraries_;
    data_t data_;
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
        tree.emplace(key_t(&slot.indexes[0], &slot.indexes[(int)slot.width]), i);
    }    
}

inline
int AlleleDepths::match_indexes_t::operator()(const key_t &rng) const {
    auto it = tree.find(rng);
    return (it != tree.end()) ? it->second : -1;
}

inline
int AlleleDepths::MatchAlleles(char *str) {
    assert(str != nullptr);
    // If string is not of length 1, return unknown
    if(str[0] == '\0' || str[1] != '\0') {
        return -1;
    }
    // code this as a switch and let the compiler optimize
    switch(str[0]) {
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'G':
    case 'g':
        return 2;
    case 'T':
    case 't':
        return 3;
    case 'N':
    case 'n':
        return 4;
    default:
        return -1;
    }
    return -1;
}

inline
uint8_t raw_counts_to_color(const depth_t& d) {
    int rank = 0;
    int zero = 0;
    for(int i=0;i<4;++i) {
        int r = 0;
        // measure the rank of the allele by comparing it to all other alleles.
        // r = 0 is the max, r = 3 is the min. Break ties by choosing the lower index.
        // make sure comparison is unsigned so we can use the high bit to signal reference
        for(int j=0;j<i;++j) {
            r += ((unsigned)d.counts[j] >= (unsigned)d.counts[i]);
        }
        for(int j=i+1;j<4;++j) {
            r += ((unsigned)d.counts[j] > (unsigned)d.counts[i]);
        }
        // rank holds the allele order in the lowest 8 bits.
        rank |= (i << (2*r));
        // zero keeps track of which alleles had a count of 0
        zero |= (d.counts[i] != 0 ? 3 : 0) << (2*r);
    }
    // xor against the maximal value, and zero out any sizes that had zero depth
    int ref = rank & 3;
    int ret = ref | (ref << 2);
    ret = ret | (ret << 4);
    ret = ((ret ^ rank) & zero);

    // restore the max reference
    ret = ret | ref;

    assert(0 <= ret && ret < 256);
    char color = AlleleDepths::hash_to_color[ret];
    assert(color != -1);

    return color;
}

struct stats_t {
    std::vector<int> node_dp;
    std::vector<int> total_depths;
    int dp;
    double log_null;
    uint8_t color;
};

inline
void calculate_stats(const RawDepths& d, int ref_index, stats_t *stats) {
    assert(stats != nullptr);

    // Measure the total depths, per-lib and overall
    int dp = 0;
    depth_t total;
    stats->node_dp.clear();
    for(auto && a : d) {
        int n = 0;
        for(int i=0;i<4;++i) {
            total.counts[i] += a.counts[i];
            n += a.counts[i];
        }
        dp += n;
        stats->node_dp.push_back(n);
    }
    stats->total_depths.assign(std::begin(total.counts), std::end(total.counts));
    stats->dp = dp;

    // calculate the log-likelihood of the null hypothesis that all reads come from binomial
    double log_null = 0.0;
    if(dp > 0) {
        log_null = -dp*log10(dp);
        for(int i=0;i<4;++i) {
            if(total.counts[i] > 0) {
                log_null += total.counts[i]*log10(total.counts[i]);
            }
        }
    }
    stats->log_null = log_null;

    int first_is_N = 0;
    if(ref_index < 4) {
        total.counts[ref_index] = utility::set_high_bit(total.counts[ref_index]);
    } else {
        first_is_N = 64;
    }
    stats->color = raw_counts_to_color(total)+first_is_N;
}

inline
void calculate_stats(const AlleleDepths& d, stats_t *stats) {
    assert(stats != nullptr);

    // Copy color
    stats->color = d.color();

    // Measure the total depths, per-lib and overall
    int dp = 0;
    std::vector<int> total(d.num_nucleotides(), 0);

    stats->node_dp.clear();
    for(size_t lib=0; lib < d.num_libraries(); ++lib) {
        int n = 0;
        for(size_t j=0; j < d.num_nucleotides(); ++j) {
            total[j] += d(lib,j);
            n += d(lib,j);
        }
        dp += n;
        stats->node_dp.push_back(n);
    }
    stats->total_depths = total;
    stats->dp = dp;

    // calculate the log-likelihood of the null hypothesis that all reads come from binomial
    double log_null = 0.0;
    if(dp > 0) {
        log_null = -dp*log10(dp);
        for(auto &&a : total) {
            if(a > 0) {
                log_null += a*log10(a);                
            }
        }
    }
    stats->log_null = log_null;
}


} // namespace pileup
} // namespace dng

#endif
