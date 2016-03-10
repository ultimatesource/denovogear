/*
 * Copyright (c) 2016 Reed A. Cartwright
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
#define BOOST_SPIRIT_USE_PHOENIX_V3 1

#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>

#include <dng/regions.h>

#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/fusion.hpp>
#include <boost/phoenix/scope/local_variable.hpp>

using namespace dng::regions;

namespace qi = boost::spirit::qi;
namespace phoenix = boost::phoenix;

BOOST_FUSION_ADAPT_STRUCT(dng::regions::detail::raw_parsed_region_t,
    (int,target),(int,beg),(int,end));

// uint_parser<unsigned, 10, 1, 3> uint3_p;        //  1..3 digits
// uint_parser<unsigned, 10, 3, 3> uint3_3_p;      //  exactly 3 digits
// test_parser("12,345,678", uint3_p >> *(',' >> uint3_3_p));


struct make_region_impl {
    typedef detail::raw_parsed_region_t result_type;

    result_type operator()(std::string s, int b, int e) const {
        return {s,b,e};
    }
    result_type operator()(std::string s, boost::optional<std::pair<int,int>> p) const {
        if(p) {
            return {s, (*p).first, (*p).second};
        }
        return {s,1,INT_MAX};
    }
};
const phoenix::function<make_region_impl> make_region;

template <typename Iterator, typename Skipper>
struct regions_grammar :
qi::grammar<Iterator, detail::raw_parsed_regions_t, Skipper> {
    typedef std::pair<int,int> pos_t;
    regions_grammar() : regions_grammar::base_type(start) {
        using qi::uint_; using qi::char_; using qi::as_string;
        using qi::attr; using qi::lexeme; using qi::no_skip;
        using qi::lit;
        using namespace qi::labels;

        start = region % no_skip[wsp];
        region = (label >> -(':' >> position))[_val = make_region(_1,_2)];
        position %= 
            ((pos[_a = _1] >> ( ('-' >> ( pos | attr(INT_MAX)))
                              | ('+' >> pos)[_1 = _a+_1-1]
                              | attr(_a)
                              ))
            |(attr(1) >> '-' >> pos)
            );
        pos = sudouble | (uint_ >> !lit(',')) | sepint;
        // We support the format of targets specified in Sam files, except that ':' are disallowed
        label = as_string[lexeme[char_("!-)+-9;<>-~") >> *char_("!-9;-~")]]; 
        sepint = lexeme[uint3_p[_val=_1] >> *(',' >> uint3_3_p[_val = 1000*_val+_1])];
    }

    qi::real_parser<double, qi::strict_ureal_policies<double>> sudouble;
    qi::uint_parser<unsigned, 10, 1, 3> uint3_p;
    qi::uint_parser<unsigned, 10, 3, 3> uint3_3_p;

    qi::rule<Iterator, detail::raw_parsed_regions_t, Skipper> start;
    qi::rule<Iterator, detail::raw_parsed_region_t, Skipper> region;
    qi::rule<Iterator, std::string(), Skipper> label;
    qi::rule<Iterator, int(), Skipper> pos, pos_end;
    qi::rule<Iterator, std::pair<int,int>(), qi::locals<int>, Skipper> position;
    qi::rule<Iterator, int(), Skipper> sepint;

    Skipper wsp;
};


std::pair<detail::raw_parsed_regions_t,bool> dng::regions::detail::parse_regions(const std::string &text) {
    namespace ss = boost::spirit::standard;
    using ss::space;
    regions_grammar<std::string::const_iterator, ss::space_type> parser_grammar;

    std::string::const_iterator first = text.begin();
    raw_parsed_regions_t regions;
    bool success = qi::phrase_parse(first, text.end(), parser_grammar, space, regions);
    success = success && (first == text.end());
    return {regions, success};
}

hts::bam::regions_t dng::regions::bam_parse(const std::string &text, const hts::bam::File &file) {
    using hts::bam::region_t;
    using hts::bam::regions_t;
    // Parse string
    auto raw_regions = dng::regions::detail::parse_regions(text);
    if(!raw_regions.second) {
        throw(std::runtime_error("Parsing of regions string failed."));
    }
    // Return quickly if there are no regions
    if(raw_regions.first.empty()) {
        return {};
    }
    // Convert raw regions to bam regions
    regions_t value;
    for(auto &&r : raw_regions.first) {
        int tid = file.TargetNameToID(r.target.c_str());
        int beg = r.beg-1;
        int end = r.end;
        if(tid < 0 ) {
            throw(std::runtime_error("Unknown contig name: '" + r.target + "'"));
        }
        if(beg < 0 || end < 0 || beg >= end) {
            throw(std::runtime_error("Invalid region coordinates: '" +
                r.target + ":" + std::to_string(beg+1) + "-" + std::to_string(end) + "'" ));
        }
        value.push_back({tid,beg,end});
    }
    // sort regions in increasing order
    std::sort(value.begin(), value.end(), [](const region_t &a, const region_t &b) {
        return std::tie(a.tid, a.beg, a.end) < std::tie(b.tid, b.beg, b.end);
    });
    // merge overlapping ranges
    auto overlapping = [](const region_t &a, const region_t &b) {
        return a.tid == b.tid && b.beg <= a.end; // assumes sorted range
    };

    auto first = value.begin();
    auto last  = value.end();
    auto result = first;
    while(++first != last) {
        if(overlapping(*result, *first)) {
            result->end = first->end;
        } else if(++result != first) {
            *result = std::move(*first);
        }
    }
    ++result;
    value.erase(result,value.end());
    return value;
}
