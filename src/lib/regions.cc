/*
 * Copyright (c) 2016 Reed A. Cartwright
 * Copyright (c) 2016 Juan Jose Garcia Mesa
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *           Juan Jose Garcia Mesa <jgarc111@asu.edu>
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
#include <vector>
#include <fstream>

#include <dng/regions.h>
#include <dng/io/utility.h>

#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/phoenix/stl/container.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/fusion.hpp>
#include <boost/phoenix/scope/local_variable.hpp>

using namespace dng::regions;

namespace qi = boost::spirit::qi;
namespace phoenix = boost::phoenix;

struct make_fragment_impl {
    typedef contig_fragment_t result_type;

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
const phoenix::function<make_fragment_impl> make_fragment;

template <typename Iterator>
struct regions_grammar :
qi::grammar<Iterator, contig_fragments_t(), boost::spirit::standard::space_type> {
    using pos_t = std::pair<int,int>;
    using skipper_type = boost::spirit::standard::space_type;
    regions_grammar() : regions_grammar::base_type(start) {
        using namespace qi;
        using namespace qi::labels;

        start = range % no_skip[wsp];
        range = (label >> -(':' >> position))[_val = make_fragment(_1,_2)];
        position %=
            ((pos[_a = _1] >> ( ('-' >> ( pos | attr(INT_MAX)))
                              | ('+' >> pos)[_1 = _a+_1-1]
                              | attr(_a)
                              ))
            |(attr(1) >> '-' >> pos)
            );
        pos = sudouble | (uint_ >> !lit(',')) | sepint;
        // We support the format of target contigs specified in Sam files, except that ':' are disallowed
        label = as_string[lexeme[char_("!-)+-9;<>-~") >> *char_("!-9;-~")]]; 
        sepint = lexeme[uint3_p[_val=_1] >> *(',' >> uint3_3_p[_val = 1000*_val+_1])];
    }

    qi::real_parser<double, qi::strict_ureal_policies<double>> sudouble;
    qi::uint_parser<unsigned, 10, 1, 3> uint3_p;
    qi::uint_parser<unsigned, 10, 3, 3> uint3_3_p;

    qi::rule<Iterator, contig_fragments_t(), skipper_type> start;
    qi::rule<Iterator, contig_fragment_t(), skipper_type> range;
    qi::rule<Iterator, std::string(), skipper_type> label;
    qi::rule<Iterator, int(), skipper_type> pos, pos_end;
    qi::rule<Iterator, std::pair<int,int>(), qi::locals<int>, skipper_type> position;
    qi::rule<Iterator, int(), skipper_type> sepint;

    skipper_type wsp;
};

boost::optional<contig_fragments_t> dng::regions::parse_contig_fragments_from_regions(const std::string &text) {
    namespace ss = boost::spirit::standard;
    using ss::space; // must match the space_type in parser_grammar
    regions_grammar<std::string::const_iterator> parser_grammar;

    std::string::const_iterator first = text.begin();
    contig_fragments_t fragments;
    bool success = qi::phrase_parse(first, text.end(), parser_grammar, space, fragments);
    success = success && (first == text.end());
    if(success) {
        return fragments;
    }
    return boost::none;
}

inline
dng::regions::ranges_t optimize_ranges(dng::regions::ranges_t value) {
    using dng::regions::ranges_t;
    using dng::regions::range_t;
     // sort regions in increasing order
    std::sort(value.begin(), value.end(), [](const range_t &a, const range_t &b) {
        return std::tie(a.beg, a.end) < std::tie(b.beg, b.end);
    });
    // merge overlapping ranges
    auto overlapping = [](const range_t &a, const range_t &b) {
        return b.beg <= a.end; // assumes sorted range
    };

    auto first = value.begin();
    auto last  = value.end();
    if(first != last) {
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
    }
    return value;
}

ranges_t dng::regions::convert_fragments_to_ranges(const contig_fragments_t& fragments, const ContigIndex& index, bool one_indexed) {
    // Convert fragments to ranges
    ranges_t value;
    int tid = 0;

    for(auto &&r : fragments) {
        tid = index.NameToId(r.contig_name, tid);
        if(tid < 0 ) {
            throw std::invalid_argument("Unknown contig name: '" + r.contig_name + "'");
        }
        int beg = r.beg-(one_indexed);
        int end = r.end;
        if( !(0 <= beg && beg <= index.contig(tid).length) ||
            !(0 <= end && end <= index.contig(tid).length) ||
            !(beg < end)
          ) {
            throw std::invalid_argument("Invalid region coordinates: '" +
                r.contig_name + ":" + std::to_string(beg+1) + "-" + std::to_string(end) + "'" );
        }
        value.push_back({tid,beg,end});
    }

    return optimize_ranges(std::move(value));
}

ranges_t dng::regions::parse_regions(const std::string &text, const ContigIndex& index) {
    // Parse string
    if(auto fragments = parse_contig_fragments_from_regions(text)) {
        return convert_fragments_to_ranges(*fragments, index, true);    
    } else {
        throw std::invalid_argument("Parsing of regions failed.");
    }    
}

template <typename Iterator>
struct bed_grammar :
qi::grammar<Iterator, contig_fragments_t(), boost::spirit::standard::blank_type> {
    using pos_t = std::pair<int,int>;
    using skipper_type = boost::spirit::standard::blank_type;
    bed_grammar() : bed_grammar::base_type(start) {
        using namespace qi;
        using namespace qi::labels;
        using boost::phoenix::push_back;

        start = (line[push_back(_val,_1)] | comment | eps) % eol;
        line = (fragment >> *(char_ - eol) );
        comment = ('#' >> *(char_ - eol));
        fragment =  (label >> uint_ >> uint_)[_val = make_fragment(_1, _2, _3)]; 
        // We support the format of target contigs specified in Sam files, except that ':' are disallowed
        label = as_string[lexeme[char_("!-)+-9;<>-~") >> *char_("!-9;-~")]];
    }

    qi::rule<Iterator, contig_fragments_t(), skipper_type> start;
    qi::rule<Iterator, contig_fragment_t(), skipper_type> line;
    qi::rule<Iterator, void(), skipper_type> comment;
    qi::rule<Iterator, contig_fragment_t(), skipper_type> fragment;
    qi::rule<Iterator, std::string(), skipper_type> label;
    qi::rule<Iterator, int(), skipper_type> pos, pos_end;
};

boost::optional<contig_fragments_t> dng::regions::parse_contig_fragments_from_bed(const std::string &text) {
    namespace ss = boost::spirit::standard;
    using ss::blank; // must match the space_type in parser_grammar
    bed_grammar<std::string::const_iterator> parser_grammar;
    std::string::const_iterator first = text.begin();
    contig_fragments_t fragments;
    bool success = qi::phrase_parse(first, text.end(), parser_grammar, blank, fragments);
    success = success && (first == text.end());
    if(success) {
        return fragments;
    }
    return boost::none;
}

dng::regions::ranges_t dng::regions::parse_bed(const std::string &text, const ContigIndex& index) {
    // Parse string
    if(auto fragments = parse_contig_fragments_from_bed(text)) {
        return convert_fragments_to_ranges(*fragments, index, false);    
    } else {
        throw std::invalid_argument("Parsing of bed failed.");
    }    
}
