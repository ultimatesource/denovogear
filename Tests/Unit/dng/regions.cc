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

#define BOOST_TEST_MODULE dng::regions

#include <dng/regions.h>

#include "../testing.h"
#include <dng/hts/hts.h>

#include <fstream>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/mem_fn.hpp>

using namespace dng::regions;
using dng::detail::make_test_range;
using hts::detail::make_data_url;

BOOST_AUTO_TEST_CASE(test_region_parsing) {
    auto test_fail = [](std::string region_string) -> void {
    BOOST_TEST_CONTEXT("region_string='" << region_string << "'") {
        auto result = parse_contig_fragments_from_regions(region_string);
        BOOST_CHECK(!result);
    }};

    test_fail("chr1:-");
    test_fail("chr1:");
    test_fail("chr1:1-1a1");

    auto test = [](std::string region_string, contig_fragments_t expected) -> void {
    BOOST_TEST_CONTEXT("region_string='" << region_string << "'") {
        using boost::adaptors::transformed;

        auto result = parse_contig_fragments_from_regions(region_string);
        BOOST_CHECK(result);

        auto test_targets = make_test_range(*result | transformed(boost::mem_fn(&contig_fragment_t::contig_name)));
        auto expected_targets = make_test_range(expected | transformed(boost::mem_fn(&contig_fragment_t::contig_name)));
        CHECK_EQUAL_RANGES(test_targets, expected_targets);

        auto test_begs = make_test_range(*result | transformed(boost::mem_fn(&contig_fragment_t::beg)));
        auto expected_begs = make_test_range(expected | transformed(boost::mem_fn(&contig_fragment_t::beg)));
        CHECK_EQUAL_RANGES(test_begs, expected_begs);

        auto test_ends = make_test_range(*result | transformed(boost::mem_fn(&contig_fragment_t::end)));
        auto expected_ends = make_test_range(expected | transformed(boost::mem_fn(&contig_fragment_t::end)));
        CHECK_EQUAL_RANGES(test_ends, expected_ends);
    }};

    test("chr1:10-1000", {{"chr1",10,1000}});
    test("chr1:10-1e3", {{"chr1",10,1000}});
    test("chr1:10-1.2e3", {{"chr1",10,1200}});
    test("chr1 : 10 - 1000", {{"chr1",10,1000}});
    test("chr1:10-",     {{"chr1",10,INT_MAX}});
    test("chr1:-1000",   {{"chr1",1,1000}});
    test("chr1:1000",    {{"chr1",1000,1000}});
    test("chr1:1,000",    {{"chr1",1000,1000}});
    test("chr1:1.0e3",    {{"chr1",1000,1000}});
    test("chr1",         {{"chr1",1,INT_MAX}});
    test("chr1:100,001 -1,001,001",    {{"chr1",100001,1001001}});

    test("chr1:10+1000", {{"chr1",10,1009}});

    test("chr1:10-1000 3:60-100", {{"chr1",10,1000},{"3",60,100}});
    test("chr1:10-1000  3", {{"chr1",10,1000},{"3",1,INT_MAX}});
    test("\nchr1\t:\n10\f-\v1000\t\n \f\v3 ", {{"chr1",10,1000},{"3",1,INT_MAX}});
}

const char bamtext[]=
"@HD\tVN:1.4\tGO:none\tSO:coordinate\n"
"@SQ\tSN:1\tLN:249250621\n"
"@SQ\tSN:2\tLN:243199373\n"
;

BOOST_AUTO_TEST_CASE(test_region_parsing_with_bam) {
    ContigIndex index;
    index.AddContig(contig_t{std::string{"1"}, 249250621});
    index.AddContig(contig_t{std::string{"2"}, 243199373});


    auto test = [&](std::string region_string, ranges_t expected) -> void {
    BOOST_TEST_CONTEXT("region_string='" << region_string << "'") {
        using boost::adaptors::transformed;

        ranges_t test;
        BOOST_REQUIRE_NO_THROW(test = parse_regions(region_string, index));

        auto test_begs = make_test_range(test | transformed(boost::mem_fn(&range_t::beg)));
        auto expected_begs = make_test_range(expected | transformed(boost::mem_fn(&range_t::beg)));
        CHECK_EQUAL_RANGES(test_begs, expected_begs);

        auto test_ends = make_test_range(test | transformed(boost::mem_fn(&range_t::end)));
        auto expected_ends = make_test_range(expected | transformed(boost::mem_fn(&range_t::end)));
        CHECK_EQUAL_RANGES(test_ends, expected_ends);
    }};

    // Check normal 
    test("1:1-1000", {{0,0,1000}});
    test("1:1-1000 2:1", {{0,0,1000},{1,0,1}});
    test("1 2", {{0,0,249250621},{1,0,243199373}});
    // Check sorting
    test("2:100-120 2:21-30 2:1-20 1:1-5 1:7-10", {{0,0,5},{0,6,10},{1,0,30},{1,99,120}});
    // Check merging
    test("1:1-1000 1:900-2000", {{0,0,2000}});
    test("2:10-200 2:100-300 1:1-1000 2:200-400", {{0,0,1000},{1,9,400}});

    // Check expected failures
    auto test_fail = [&](std::string region_string) -> void {
    BOOST_TEST_CONTEXT("region_string='" << region_string << "'") {
        ranges_t test;
        BOOST_CHECK_THROW(test = parse_regions(region_string, index), std::invalid_argument);
    }};

    test_fail("xyz");
    test_fail("1:1000-900");
    test_fail("1:1000-999");
    // Check empty file
    test_fail("");
}

BOOST_AUTO_TEST_CASE(test_bed_parsing) {
    ContigIndex index;
    index.AddContig(contig_t{std::string{"1"}, 249250621});
    index.AddContig(contig_t{std::string{"2"}, 243199373});
    
    auto test = [&](std::string bed_string, ranges_t expected) -> void {
    BOOST_TEST_CONTEXT("bed_string='" << bed_string << "'") {
        using boost::adaptors::transformed;

        ranges_t test;
        BOOST_REQUIRE_NO_THROW(test = parse_bed(bed_string, index));

        auto test_begs = make_test_range(test | transformed(boost::mem_fn(&range_t::beg)));
        auto expected_begs = make_test_range(expected | transformed(boost::mem_fn(&range_t::beg)));
        CHECK_EQUAL_RANGES(test_begs, expected_begs);

        auto test_ends = make_test_range(test | transformed(boost::mem_fn(&range_t::end)));
        auto expected_ends = make_test_range(expected | transformed(boost::mem_fn(&range_t::end)));
        CHECK_EQUAL_RANGES(test_ends, expected_ends);
    }};

    // Check normal
    test("1\t0\t1000\n", {{0,0,1000}});
    test("1\t0\t1000\n2\t0\t1", {{0,0,1000},{1,0,1}});
    // Check sorting
    test("2\t99\t120\n2\t20\t30\n2\t0\t20\n1\t0\t5\n1\t6\t10", {{0,0,5},{0,6,10},{1,0,30},{1,99,120}});
    // Check merging
    test("1\t0\t1000\n1\t899\t2000", {{0,0,2000}});
    test("2\t9\t200\n2\t99\t300\n1\t0\t1000\n2\t199\t400", {{0,0,1000},{1,9,400}});
    // Check Comments and blank lines
    test("#comment\n\n\n#comment\ttest\n#comment\ttest\tanothertest\n1\t0\t1000\n", {{0,0,1000}});
    // Check empty file
    test("", {});
}
