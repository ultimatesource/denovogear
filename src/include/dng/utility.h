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

#include <boost/spirit/include/support_ascii.hpp>
#include <boost/spirit/include/qi_real.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/qi_list.hpp>
#include <boost/spirit/include/qi_char.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/algorithm/string/predicate.hpp>

#include <dng/detail/unit_test.h>

namespace dng {
namespace utility {

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

// TODO: make the separator more generic
template<typename S>
std::pair<std::vector<double>, bool> parse_double_list(const S &str,
        char sep = ',', size_t sz_hint = 4) {
    namespace qi = boost::spirit::qi;
    std::vector<double> f;
    f.reserve(sz_hint);
    boost::spirit::standard::space_type space;
    auto b = boost::begin(str);
    auto e = boost::end(str);
    bool r = qi::phrase_parse(b, e, qi::double_ % sep, space, f);
    return {f, (r &&b == e) };
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
inline std::pair<std::string, std::string> extract_file_type(const std::string &path) {
    if(path.empty())
        return {};
    std::locale loc;

    auto last = path.length();
    decltype(last) first = 0;
    while(first < path.length() && std::isspace(path[first],loc)) {
        ++first;
    }
    if(first == last) {
        return {};
    }
    while(std::isspace(path[last-1],loc)) {
        --last;
    }

    auto x = last-1;
    for(auto u = first; u < last; ++u) {
        if(path[u] == ':' && u > first+1) { // u > 1 skips windows drive letters
            return {path.substr(first, u-first), path.substr(u+1,last-(u+1))};
        }
        if(path[u] == '.' && u > first) { // u > 0 skips unix .hidden files
            x = u;
        }
    }
    return {path.substr(x + 1, last - (x + 1)), path.substr(first,last-first)};
}

}

using utility::location_t;

} // namespace dng::utility

#endif
