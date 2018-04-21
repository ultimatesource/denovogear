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
#include <limits>
#include <array>
#include <type_traits>

#include <boost/spirit/include/support_ascii.hpp>
#include <boost/spirit/include/qi_real.hpp>
#include <boost/spirit/include/qi_int.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/qi_list.hpp>
#include <boost/spirit/include/qi_char.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/tokenizer.hpp>

#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>

#include <boost/range/iterator.hpp>
#include <boost/range/distance.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/as_literal.hpp>
#include <boost/range/reference.hpp>
#include <boost/range/size.hpp>

#include <dng/detail/unit_test.h>

namespace dng {
namespace utility {

using StringSet = boost::container::flat_set<std::string>;
using StringMap = boost::container::flat_map<std::string, size_t>;

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

template<typename Range, typename Value>
inline
size_t find_position(const Range& r, const Value& v) {
    return boost::distance(boost::find<boost::return_begin_found>(r, v));
}

template<typename Range, typename Pred>
inline
size_t find_position_if(const Range& r, Pred pred) {
    return boost::distance(boost::find_if<boost::return_begin_found>(r, pred));
}

template<typename A, typename B>
std::size_t key_switch(const A &ss, const B &keys) {
    std::size_t ret = find_position_if(keys,
        [&](typename boost::range_reference<const B>::type k) {
            return boost::algorithm::istarts_with(k, ss);
        });
    return (ret != boost::size(keys)) ? ret : static_cast<std::size_t>(-1);
}

template<typename A, typename B>
std::size_t key_switch_iequals(const A &ss, const B &keys) {
    std::size_t ret = find_position_if(keys,
        [&](typename boost::range_reference<const B>::type k) {
            return boost::algorithm::iequals(k, ss);
        });
    return (ret != boost::size(keys)) ? ret : static_cast<std::size_t>(-1);
}

template<typename A, typename B>
typename boost::range_reference<const B>::type
key_switch_tuple(const A &ss, const B &keys,
    typename boost::range_reference<const B>::type &default_value) {
    auto it = boost::find_if<boost::return_found>(keys,
        [&](typename boost::range_reference<const B>::type k) {
            return boost::algorithm::istarts_with(std::get<0>(k), ss);
        });
    return (it != boost::end(keys)) ? *it : default_value;
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
        throw std::invalid_argument("Unable to parse nuc-freq option. "
                                 "It must be a comma separated list of floating-point numbers.");
    }
    if(f.first.size() != 4) {
        throw std::invalid_argument("Wrong number of values passed to nuc-freq. "
                                 "Expected 4; found " + std::to_string(f.first.size()) + ".");
    }
    std::array<double,4> freqs;
    std::copy(f.first.begin(), f.first.end(), &freqs[0]);

    // Scale frequencies and makes sure they sum to 1.
    double freqs_sum = 0.0;
    for(int i=0;i<freqs.size();++i) {
        if(freqs[i] < 0.0) {
            throw std::invalid_argument("Nucleotide frequencies must be non-negative. "
                "Parsing of '" + str + "' generated a frequency of " + std::to_string(freqs[i])
                + " in position " + std::to_string(i) + "."
                );
        }
        freqs_sum += freqs[i];
    }
    for(auto && a : freqs) {
        a = a/freqs_sum;
    }

    return freqs;
}

template<typename T>
std::string to_pretty(const T &value) {
    namespace karma = boost::spirit::karma;
    std::string str;
    str.reserve(16);
    if(karma::generate(std::back_inserter(str), value)) {
        return str;
    }
    return {};
}

/*
  Phred scaled numbers: -10.0*log10(a)
  
  DBL_MIN -> 3076.52655568588761525
  DBL_SUBNORM_MIN -> 3233.06215343115809446
  1-(1-2^-53) -> 159.54589770191000753
*/

// -10*log10(a)
inline double phred(double a) {
    // Use an if to prevent -0.0
    return (a == 1.0) ? 0.0 : -10.0 * std::log10(a);
}

// -10*log10(1+p)
inline double phred1p(double p) {
    // Use an if to prevent -0.0
    return (p == 0.0) ? 0.0 : (-10.0/M_LN10) * std::log1p(p);
}

template<typename T>
inline T lphred(double a, T m = std::numeric_limits<T>::max()) {
    double q = std::round( -10.0 * std::log10(a) );
    return (q > m) ? m : static_cast<T>(q);
}

template<typename T>
inline T lphred1p(double p, T m = std::numeric_limits<T>::max()) {
    double q = std::round( (-10.0/M_LN10) * std::log1p(p) );
    return (q > m) ? m : static_cast<T>(q);
}

// extracts extension and filename from both file.foo and ext:file.foo
// returns {ext, file.foo}
// trims whitespace as well
std::pair<std::string, std::string> extract_file_type(const char *path);

inline
std::pair<std::string, std::string> extract_file_type(const std::string &path) {
    return extract_file_type(path.c_str());
}

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
// converts an input file to category and will throw if value is not supported
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

namespace detail {
    using token_function = boost::char_separator<char>;
    template<typename Range>
    using char_tokenizer = boost::tokenizer<token_function, typename boost::range_iterator<const Range>::type>;
};

template<typename Range>
inline
detail::char_tokenizer<Range>
make_tokenizer(const Range& text, const char *sep = "\t", const char *eol = "\n") {
    detail::token_function f(sep, eol, boost::keep_empty_tokens);
    return {boost::as_literal(text),f};
}

template<typename Range>
inline
detail::char_tokenizer<Range>
make_tokenizer_dropempty(const Range& text, const char *sep = "\t", const char *eol = "\n") {
    detail::token_function f(sep, eol, boost::drop_empty_tokens);
    return {boost::as_literal(text),f};
}

template<typename T>
inline T set_high_bit(T x) {
  T mask = 1;
  mask <<= std::numeric_limits<typename std::make_signed<T>::type>::digits;
  return x | mask;
}

// this make_array implementation taken from
// https://gist.github.com/klmr/2775736
template <typename... T>
inline
constexpr auto make_array(T&&... values) ->
        std::array<
            typename std::decay<
                typename std::common_type<T...>::type>::type,
            sizeof...(T)> {
    return std::array<
        typename std::decay<
            typename std::common_type<T...>::type>::type,
        sizeof...(T)>{std::forward<T>(values)...};
}

} // namespace dng::utility

using utility::location_t;

} // namespace dng

#endif
