/*
 * Copyright (c) 2014-2016 Reed A. Cartwright
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
#ifndef DNG_UTILITY_H
#define DNG_UTILITY_H

#include <tuple>
#include <cmath>
#include <locale>
#include <cstdint>
#include <climits>
#include <chrono>


#include <boost/spirit/include/support_ascii.hpp>
#include <boost/spirit/include/qi_real.hpp>
#include <boost/spirit/include/qi_int.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/qi_list.hpp>
#include <boost/spirit/include/qi_char.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/algorithm/string/predicate.hpp>

#include <dng/detail/unit_test.h>

namespace dng {
namespace utility {

template<typename E, typename = typename std::enable_if<std::is_enum<E>::value>::type>
struct EnumFlags {
    typedef E enum_t;
    typedef typename std::underlying_type<E>::type int_t;
    static_assert(std::is_enum<enum_t>::value, "EnumFlags can only wrap an enum.");

    int_t value;
    
    EnumFlags(int_t v) : value{v} { }
    EnumFlags(enum_t v) : value{static_cast<int_t>(v)} { }

    EnumFlags& operator=(enum_t v) {
        value = static_cast<int_t>(v);
        return *this;
    }

    operator int_t() { return value; }

    EnumFlags operator|(enum_t v) {
         return {value | static_cast<int_t>(v)};
    }

    EnumFlags& operator|=(enum_t v) {
         value |= static_cast<int_t>(v);
         return *this;
    }

    EnumFlags operator&(enum_t v) {
        return {value & static_cast<int_t>(v)};
    }
};

template<typename E>
EnumFlags<E> operator|(E a, E b) {
    return EnumFlags<E>{a} | b;
}

typedef int64_t location_t;

constexpr location_t LOCATION_MAX = ((static_cast<location_t>(INT32_MAX) << 32) | INT32_MAX);

inline location_t make_location(int t, int p) {
    assert((t & INT32_MAX) == t && (p & INT32_MAX) == p);
    return (static_cast<location_t>(t) << 32) | p;
}
inline int location_to_contig(location_t u) {
    int x = static_cast<int>(u >> 32);
    assert((x & INT32_MAX) == x);
    return x;
}
inline int location_to_position(location_t u) {
    int x = static_cast<int>(u & UINT32_MAX); // yes, UINT32_MAX
    assert((x & INT32_MAX) == x);
    return x;
}

template<class A, class B, std::size_t N>
std::size_t key_switch(A &ss, const B(&key)[N]) {
    using boost::algorithm::istarts_with;
    for(std::size_t i = 0; i < N; ++i) {
        if(istarts_with(key[i], ss)) {
            return i;
        }
    }
    return static_cast<std::size_t>(-1);
}

template<class A, class B, std::size_t N>
const B &key_switch_tuple(A &ss, const B(&key)[N], const B &default_value) {
    using boost::algorithm::istarts_with;
    for(std::size_t i = 0; i < N; ++i) {
        if(istarts_with(std::get<0>(key[i]), ss)) {
            return key[i];
        }
    }
    return default_value;
}

template<class A, class B, std::size_t N>
std::size_t key_switch_iequals(A &ss, const B(&key)[N]) {
    using boost::algorithm::iequals;
    for(std::size_t i = 0; i < N; ++i) {
        if(iequals(key[i], ss)) {
            return i;
        }
    }
    return static_cast<std::size_t>(-1);
}

// TODO: make the separator more generic
template<typename S>
std::pair<std::vector<double>, bool> parse_double_list(const S &str,
        char sep = ',', size_t sz_hint = 4) {
    namespace qi = boost::spirit::qi;
    namespace ss = boost::spirit::standard;
    using ss::blank;

    std::vector<double> f;
    f.reserve(sz_hint);
    auto b = boost::begin(str);
    auto e = boost::end(str);
    bool r = qi::phrase_parse(b, e, qi::double_ % sep, blank, f);
    return {f, (r && b == e) };
}

template<typename S>
std::pair<std::vector<int>, bool> parse_int_list(const S &str,
        char sep = ',', size_t sz_hint = 4) {
    namespace qi = boost::spirit::qi;
    namespace ss = boost::spirit::standard;
    using ss::blank;

    std::vector<int> f;
    f.reserve(sz_hint);
    auto b = boost::begin(str);
    auto e = boost::end(str);
    bool r = qi::phrase_parse(b, e, qi::int_ % sep, blank, f);
    return {f, (r && b == e) };
}

// Parse Nucleotide Frequencies
inline
std::array<double, 4> parse_nuc_freqs(const std::string &str) {
    auto f = utility::parse_double_list(str, ',', 4);
    if(!f.second) {
        throw std::runtime_error("Unable to parse nuc-freq option. "
                                 "It must be a comma separated list of floating-point numbers.");
    }
    if(f.first.size() != 4) {
        throw std::runtime_error("Wrong number of values passed to nuc-freq. "
                                 "Expected 4; found " + std::to_string(f.first.size()) + ".");
    }
    std::array<double,4> freqs;
    std::copy(f.first.begin(), f.first.begin()+4, &freqs[0]);
    return freqs;
}

template<typename T>
std::string to_pretty(const T &value) {
    namespace karma = boost::spirit::karma;
    std::string str;
    if(karma::generate(std::back_inserter(str), value)) {
        return str;
    }
    return {};
}

inline double phred(double p) {
    // Use an if to prevent -0.0
    return (p == 1.0) ? 0.0 : -10.0 * std::log10(p);
}

template<typename T>
inline T lphred(double p, T m = std::numeric_limits<T>::max()) {
    double q = std::round(phred(p));
    return (q > m) ? m : static_cast<T>(q);
}

// extracts extension and filename from both file.ext and ext:file.foo
// returns {ext, filename.ext}
// trims whitespace as well
std::pair<std::string, std::string> extract_file_type(const std::string &path);

// a strongly-typed enum for file category
enum class FileCat {
    Unknown  = 0,
    Sequence = 1,
    Variant  = 2,
    Pileup   = 4
};
typedef EnumFlags<FileCat> FileCatSet;

// converts an extension to a file category
FileCat file_category(const std::string &ext);
// converts and input file to category and will throw if value is not supported
FileCat input_category(const std::string &in, FileCatSet mask, FileCat def = FileCat::Unknown);

// create a timestamp that contains the date and epoch
std::pair<std::string, std::string> timestamp();

// create a timestamp in "Date=xxx,Epoch=xxx" format
std::string vcf_timestamp();

template<typename V, typename A>
inline
std::string vcf_command_line_text(const char *arg,
                                  const std::vector<V, A> &val) {
    std::string str;
    for(auto && a : val) {
        str += std::string("--") + arg + '=' + dng::utility::to_pretty(a) + ' ';
    }
    str.pop_back();
    return str;
}


template<typename VAL>
inline
std::string vcf_command_line_text(const char *arg, VAL val) {
    return std::string("--") + arg + '=' + dng::utility::to_pretty(val);
}

template<>
inline
std::string vcf_command_line_text(const char *arg, std::string val) {
    return std::string("--") + arg + "=\'" + val + "\'";
}

} // namespace dng::utility

using utility::location_t;

} // namespace dng

#endif
