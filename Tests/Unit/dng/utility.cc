/*
 * Copyright (c) 2016-2017 Reed A. Cartwright
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

#define BOOST_TEST_MODULE dng::utility

#include <dng/utility.h>

#include "../testing.h"

#include <vector>
#include <string>
#include <utility>

#include <boost/range/algorithm/copy.hpp>
#include <iterator>

namespace dng { namespace utility {
    // pull in custom io operator for classes in this namespace
    using dng::detail::io::operator<<;
}}

using namespace std;
using namespace dng::utility;
using namespace dng::detail::io;
using dng::detail::make_test_range;

BOOST_AUTO_TEST_CASE(test_locations) {
    BOOST_CHECK_EQUAL(LOCATION_MAX, 0x7FFFFFFF7FFFFFFFll);

    auto test = [](int t, int p) -> void {
    BOOST_TEST_CONTEXT("t=" << t << ", p=" << p) {
        BOOST_REQUIRE_GE(t, 0);
        BOOST_REQUIRE_GE(p, 0);
        BOOST_REQUIRE_LE(t, INT32_MAX);
        BOOST_REQUIRE_LE(p, INT32_MAX);

        location_t test_location = make_location(t,p);
        int64_t expected_location = ((int64_t)t << 32) | p;
        BOOST_CHECK_EQUAL(test_location, expected_location);
        BOOST_CHECK_EQUAL(location_to_contig(test_location), t);
        BOOST_CHECK_EQUAL(location_to_position(test_location), p);
    }};

    test(0,0);
    test(0,1);
    test(10,2147483647);
    test(2147483647,2147483647);
}

BOOST_AUTO_TEST_CASE(test_key_switch) {
    string keys[] = {
        "0","1","2","male","female","unknown","a1","a2"
    };

    BOOST_CHECK_EQUAL(key_switch("x", keys), -1);
    BOOST_CHECK_EQUAL(key_switch("", keys), 0);
    BOOST_CHECK_EQUAL(key_switch("0", keys), 0);
    BOOST_CHECK_EQUAL(key_switch("1", keys), 1);
    BOOST_CHECK_EQUAL(key_switch("2", keys), 2);
    BOOST_CHECK_EQUAL(key_switch("m", keys), 3);
    BOOST_CHECK_EQUAL(key_switch("m", keys), 3);
    BOOST_CHECK_EQUAL(key_switch("f", keys), 4);
    BOOST_CHECK_EQUAL(key_switch("u", keys), 5);

    BOOST_CHECK_EQUAL(key_switch("a", keys), 6);
    BOOST_CHECK_EQUAL(key_switch("a1", keys), 6);
    BOOST_CHECK_EQUAL(key_switch("a2", keys), 7);
    BOOST_CHECK_EQUAL(key_switch("A1", keys), 6);
    BOOST_CHECK_EQUAL(key_switch("A2", keys), 7);

    BOOST_CHECK_EQUAL(key_switch("Male", keys), 3);
    BOOST_CHECK_EQUAL(key_switch("MALE", keys), 3);
    BOOST_CHECK_EQUAL(key_switch("MaLe", keys), 3);
    BOOST_CHECK_EQUAL(key_switch("mal", keys), 3);
    BOOST_CHECK_EQUAL(key_switch("ma", keys), 3);
    BOOST_CHECK_EQUAL(key_switch("male ", keys), -1);
}

BOOST_AUTO_TEST_CASE(test_key_switch_iequals) {
    string keys[] = {
        "0","1","2","male","female","unknown","a1","a2"
    };

    BOOST_CHECK_EQUAL(key_switch_iequals("0", keys), 0);
    BOOST_CHECK_EQUAL(key_switch_iequals("1", keys), 1);
    BOOST_CHECK_EQUAL(key_switch_iequals("2", keys), 2);
    BOOST_CHECK_EQUAL(key_switch_iequals("male", keys), 3);
    BOOST_CHECK_EQUAL(key_switch_iequals("female", keys), 4);
    BOOST_CHECK_EQUAL(key_switch_iequals("unknown", keys), 5);
    BOOST_CHECK_EQUAL(key_switch_iequals("a1", keys), 6);
    BOOST_CHECK_EQUAL(key_switch_iequals("a2", keys), 7);

    BOOST_CHECK_EQUAL(key_switch_iequals("x", keys), -1);
    BOOST_CHECK_EQUAL(key_switch_iequals("", keys), -1);
    BOOST_CHECK_EQUAL(key_switch_iequals("m", keys), -1);
    BOOST_CHECK_EQUAL(key_switch_iequals("m", keys), -1);
    BOOST_CHECK_EQUAL(key_switch_iequals("f", keys), -1);
    BOOST_CHECK_EQUAL(key_switch_iequals("u", keys), -1);

    BOOST_CHECK_EQUAL(key_switch_iequals("a", keys), -1);
    BOOST_CHECK_EQUAL(key_switch_iequals("A1", keys), 6);
    BOOST_CHECK_EQUAL(key_switch_iequals("A2", keys), 7);

    BOOST_CHECK_EQUAL(key_switch_iequals("Male", keys), 3);
    BOOST_CHECK_EQUAL(key_switch_iequals("MALE", keys), 3);
    BOOST_CHECK_EQUAL(key_switch_iequals("MaLe", keys), 3);
    BOOST_CHECK_EQUAL(key_switch_iequals("mal", keys), -1);
    BOOST_CHECK_EQUAL(key_switch_iequals("ma", keys), -1);
}

BOOST_AUTO_TEST_CASE(test_key_switch_tuple) {
    pair<string,int> keys[] = {
        {"0", 0},
        {"1", 1},
        {"2", 2},
        {"male", 1},
        {"female", 2},
        {"unknown", 0}
    };

    BOOST_CHECK_EQUAL(key_switch_tuple("", keys, keys[0]).second, 0);
    BOOST_CHECK_EQUAL(key_switch_tuple("0", keys, keys[0]).second, 0);
    BOOST_CHECK_EQUAL(key_switch_tuple("1", keys, keys[0]).second, 1);
    BOOST_CHECK_EQUAL(key_switch_tuple("2", keys, keys[0]).second, 2);
    BOOST_CHECK_EQUAL(key_switch_tuple("male", keys, keys[0]).second, 1);
    BOOST_CHECK_EQUAL(key_switch_tuple("m", keys, keys[0]).second, 1);
    BOOST_CHECK_EQUAL(key_switch_tuple("female", keys, keys[0]).second, 2);
    BOOST_CHECK_EQUAL(key_switch_tuple("f", keys, keys[0]).second, 2);
    BOOST_CHECK_EQUAL(key_switch_tuple("unknown", keys, keys[0]).second, 0);
    BOOST_CHECK_EQUAL(key_switch_tuple("MALE", keys, keys[0]).second, 1);
    BOOST_CHECK_EQUAL(key_switch_tuple("M", keys, keys[0]).second, 1);
    BOOST_CHECK_EQUAL(key_switch_tuple("Female", keys, keys[0]).second, 2);
    BOOST_CHECK_EQUAL(key_switch_tuple("F", keys, keys[0]).second, 2);
}

BOOST_AUTO_TEST_CASE(test_parse_double_list) {
    string nfreq = " 0.3, 0.2 ,0.2 , 0.3 ";
    vector<double> nfreq_a = {0.3, 0.2, 0.2, 0.3};
    char sep = ',';
    size_t size = 4;

    pair<vector<double>, bool> dup = parse_double_list(nfreq, sep, size);
    BOOST_CHECK_EQUAL(dup.second, true);
    CHECK_EQUAL_RANGES(dup.first, nfreq_a);
}

BOOST_AUTO_TEST_CASE(test_parse_int_list) {
    string nfreq = " 1, 2 ,3 , 4  ";
    vector<double> nfreq_a = {1,2,3,4};
    char sep = ',';
    size_t size = 4;

    pair<vector<int>, bool> dup = parse_int_list(nfreq, sep, size);
    BOOST_CHECK_EQUAL(dup.second, true);
    CHECK_EQUAL_RANGES(dup.first, nfreq_a);
}

BOOST_AUTO_TEST_CASE(test_to_pretty) {
    BOOST_CHECK_EQUAL(to_pretty(1), "1");
    BOOST_CHECK_EQUAL(to_pretty(1.234), "1.234");
    BOOST_CHECK_EQUAL(to_pretty("string"), "string");
    BOOST_CHECK_EQUAL(to_pretty('c'), "c");
}

BOOST_AUTO_TEST_CASE(test_extract_file_type) {
    using pair = std::pair<std::string,std::string>;

    BOOST_CHECK_EQUAL(extract_file_type("foo.bar"), pair("bar", "foo.bar"));
    BOOST_CHECK_EQUAL(extract_file_type("my:foo.bar"), pair("my", "foo.bar"));
    BOOST_CHECK_EQUAL(extract_file_type(".bar"), pair("", ".bar"));
    BOOST_CHECK_EQUAL(extract_file_type("my:.foo.bar"), pair("my",".foo.bar"));
    BOOST_CHECK_EQUAL(extract_file_type(".foo.bar"), pair("bar",".foo.bar"));
    BOOST_CHECK_EQUAL(extract_file_type(""), pair("", ""));
    BOOST_CHECK_EQUAL(extract_file_type(std::string{}), pair({},{}));
    //BOOST_CHECK(extract_file_type("C:\\foo.bar"), pair("bar"},string{"C:\\foo.bar"}));
    //BOOST_CHECK(extract_file_type("CC:C:\\foo.bar"), pair("CC"},string{"C:\\foo.bar"}));

    BOOST_CHECK_EQUAL(extract_file_type("bcf:test"), pair("bcf", "test"));
    BOOST_CHECK_EQUAL(extract_file_type("vcf:-"), pair("vcf", "-"));
    BOOST_CHECK_EQUAL(extract_file_type("vcf:"), pair("vcf", ""));
    BOOST_CHECK_EQUAL(extract_file_type("test.bcf"), pair("bcf", "test.bcf"));
    BOOST_CHECK_EQUAL(extract_file_type("test.vcf"), pair("vcf", "test.vcf"));

    BOOST_CHECK_EQUAL(extract_file_type(" \f\n\r\t\vfoo.bar \f\n\r\t\v"), pair("bar", "foo.bar"));
    BOOST_CHECK_EQUAL(extract_file_type(" \f\n\r\t\vmy:foo.bar \f\n\r\t\v"), pair("my", "foo.bar"));
    BOOST_CHECK_EQUAL(extract_file_type(" \f\n\r\t\v.bar \f\n\r\t\v"), pair("", ".bar"));
    BOOST_CHECK_EQUAL(extract_file_type(" \f\n\r\t\v"), pair({},{}));
}

BOOST_AUTO_TEST_CASE(test_file_category) {
    using namespace std;
    using dng::utility::FileCat;

    BOOST_CHECK_EQUAL(file_category(""),     FileCat::Unknown);
    BOOST_CHECK_EQUAL(file_category("txt"),  FileCat::Unknown);
    BOOST_CHECK_EQUAL(file_category("bam"),  FileCat::Sequence);
    BOOST_CHECK_EQUAL(file_category("cram"), FileCat::Sequence);
    BOOST_CHECK_EQUAL(file_category("sam"),  FileCat::Sequence);
    BOOST_CHECK_EQUAL(file_category("vcf"),  FileCat::Variant);
    BOOST_CHECK_EQUAL(file_category("bcf"),  FileCat::Variant);
    BOOST_CHECK_EQUAL(file_category("ad"),   FileCat::Pileup);
    BOOST_CHECK_EQUAL(file_category("tad"),  FileCat::Pileup);

    BOOST_CHECK_EQUAL(file_category("BAM"),  FileCat::Sequence);
    BOOST_CHECK_EQUAL(file_category("CRAM"), FileCat::Sequence);
    BOOST_CHECK_EQUAL(file_category("SAM"),  FileCat::Sequence);
    BOOST_CHECK_EQUAL(file_category("VCF"),  FileCat::Variant);
    BOOST_CHECK_EQUAL(file_category("BCF"),  FileCat::Variant);
    BOOST_CHECK_EQUAL(file_category("AD"),   FileCat::Pileup);
    BOOST_CHECK_EQUAL(file_category("TAD"),  FileCat::Pileup);
}

BOOST_AUTO_TEST_CASE(test_tokenizer) {
    const char text[] = "1\t2\t3\t4\t\n5\n\n6\t7\t8\t9\n";
    auto tokenizer = make_tokenizer(text);
    vector<string> test_tokens;
    boost::copy(tokenizer, std::back_inserter(test_tokens));
    vector<string> expected_tokens = {
        "1", "2", "3", "4", "", "\n", "5", "\n", "", "\n",
        "6", "7", "8", "9", "\n", ""
    };
    CHECK_EQUAL_RANGES(test_tokens, expected_tokens);

    auto tokenizer_de = make_tokenizer_dropempty(text);
    vector<string> test_tokens_dropempty;
    boost::copy(tokenizer_de, std::back_inserter(test_tokens_dropempty));
    vector<string> expected_tokens_dropempty = {
        "1", "2", "3", "4", "\n", "5", "\n", "\n",
        "6", "7", "8", "9", "\n"
    };
    CHECK_EQUAL_RANGES(test_tokens_dropempty, expected_tokens_dropempty);
}

BOOST_AUTO_TEST_CASE(test_set_high_bit) {
    BOOST_CHECK_EQUAL(set_high_bit<uint8_t>(0),  0x80u);
    BOOST_CHECK_EQUAL(set_high_bit<uint16_t>(0), 0x8000u);
    BOOST_CHECK_EQUAL(set_high_bit<uint32_t>(0), 0x80000000u);
    BOOST_CHECK_EQUAL(set_high_bit<uint64_t>(0), 0x8000000000000000llu);
}
