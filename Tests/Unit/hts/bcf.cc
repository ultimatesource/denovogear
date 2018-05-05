/*
 * Copyright (c) 2018 Reed A. Cartwright
 * Copyright (c) 2018 Juan J. Garcia Mesa
 * Authors: Reed A. Cartwright <reed@cartwrig.ht>
 *	    Juan J. Garcia Mesa <jgarc111@asu.edu>
 *
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distrubited in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY of FITNESS
 * FOR ANY PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE hts::bcf

#include <dng/hts/bcf.h>
#include <dng/hts/hts.h>

#include "../testing.h"

const char vcfformatint[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1,2,3,4\n"
;

const char vcfformatfloat[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Float>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1.1,2.1,3.1,4.1\n"
;

const char vcftypechar[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Character>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1,2,3,4\n"
;

const char vcftypestring[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=String>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1,2,3,4\n"
;

const char vcfupdateint[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20\n"
;

const char vcfupdatefloat[] =
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=PL,Number=G,Type=Float>\n"
    "##contig=<ID=1,length=100>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\t.\tPL\t1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0\n"
;

const char vcfinfo[] =
    "##fileformat=VCFv4.2\n"
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n"
    "##contig=<ID=1,length=100>\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"# high-quality bases\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
    "1\t1\t.\tG\tA\t.\t.\tDP=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20\t.\t.\n"
;

using namespace hts;
using namespace hts::bcf;
using hts::detail::make_data_url;
using dng::detail::make_test_range;

BOOST_AUTO_TEST_CASE(test_vcf_read_int_PL) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, const std::vector<int> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        // Set header and rec for reading PL values
        auto file = hts_open(vcfstr.c_str(), "r");
        auto header = bcf_hdr_read(file);
        auto record = bcf_init();
        bcf_read(file, header, record);
        bcf_unpack(record, BCF_UN_STR);

        int n_pl = 8;
        buffer_t<int> buffer = hts::bcf::make_buffer<int>(n_pl);

        // Get PL values
        int n = hts::bcf::get_format_numeric(header, record, "PL", &buffer, &n_pl);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    test(make_data_url(vcfformatint), {1,2,3,4});
    test(make_data_url(vcfformatfloat), {1,2,3,4});
    test(make_data_url(vcfupdateint), {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20});
    test(make_data_url(vcfupdatefloat), {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20});
}

BOOST_AUTO_TEST_CASE(test_vcf_read_float_PL) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, const std::vector<float> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        // Set header and rec for reading PL values
        auto file = hts_open(vcfstr.c_str(), "r");
        auto header = bcf_hdr_read(file);
        auto record = bcf_init();
        bcf_read(file, header, record);
        bcf_unpack(record, BCF_UN_STR);

        int n_pl = 8;
        buffer_t<float> buffer = hts::bcf::make_buffer<float>(n_pl);

        // Get PL values
        int n = hts::bcf::get_format_numeric(header, record, "PL", &buffer, &n_pl);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    test(make_data_url(vcfformatint), {1.0,2.0,3.0,4.0});
    test(make_data_url(vcfformatfloat), {1.1,2.1,3.1,4.1});
    test(make_data_url(vcfupdateint), {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0});
    test(make_data_url(vcfupdatefloat), {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0});
}

BOOST_AUTO_TEST_CASE(test_vcf_read_PL_tag_fail) {
    auto test = [](const std::string &vcfstr, const char *tag, int expected) -> void {
    BOOST_TEST_CONTEXT("tag=" << tag) {
	// Set header and rec for reading PL values
        auto file = hts_open(vcfstr.c_str(), "r");
        auto header = bcf_hdr_read(file);
        auto record = bcf_init();
        bcf_read(file, header, record);
        bcf_unpack(record, BCF_UN_STR);

        int n_pl = 8;//10*record->n_sample;
        buffer_t<float> buffer_int = hts::bcf::make_buffer<float>(n_pl);
	buffer_t<int> buffer_float = hts::bcf::make_buffer<int>(n_pl);
        // Get PL values
        int n = hts::bcf::get_format_numeric(header, record, tag, &buffer_int, &n_pl);
	BOOST_REQUIRE_EQUAL(n,expected);
	n = hts::bcf::get_format_numeric(header, record, tag, &buffer_float, &n_pl);
	BOOST_REQUIRE_EQUAL(n,expected);
    }};

    test(make_data_url(vcfformatint), "Pl",-1);
    test(make_data_url(vcfformatint), "pL",-1);
    test(make_data_url(vcfformatint), "pl",-1);
    test(make_data_url(vcfformatint), "Pl",-1);
}

BOOST_AUTO_TEST_CASE(test_vcf_read_PL_type_fail) {
    auto test = [](const std::string &vcfstr, int expected) -> void {
    BOOST_TEST_CONTEXT(vcfstr) {
	// Set header and rec for reading PL values
        auto file = hts_open(vcfstr.c_str(), "r");
        auto header = bcf_hdr_read(file);
        auto record = bcf_init();
        bcf_read(file, header, record);
        bcf_unpack(record, BCF_UN_STR);

        int n_pl = 10*record->n_sample;
        buffer_t<float> buffer_int = hts::bcf::make_buffer<float>(n_pl);
	buffer_t<int> buffer_float = hts::bcf::make_buffer<int>(n_pl);
        // Get PL values
        int n = hts::bcf::get_format_numeric(header, record, "PL", &buffer_int, &n_pl);
	BOOST_REQUIRE_EQUAL(n,expected);
	n = hts::bcf::get_format_numeric(header, record, "PL", &buffer_float, &n_pl);
	BOOST_REQUIRE_EQUAL(n,expected);
    }};

    test(make_data_url(vcftypechar), -2);
    test(make_data_url(vcftypestring), -2);
}

BOOST_AUTO_TEST_CASE(test_vcf_info_DP) {
    int counter = 0;
    auto test = [&](const std::string &vcfstr, const std::vector<int> &expected) -> void {
    BOOST_TEST_CONTEXT("counter=" << ++counter) {
        // Set header and rec for reading PL values
        auto file = hts_open(vcfstr.c_str(), "r");
        auto header = bcf_hdr_read(file);
        auto record = bcf_init();
        bcf_read(file, header, record);
        bcf_unpack(record, BCF_UN_STR);

        int n_pl = 8;
        buffer_t<int> buffer = hts::bcf::make_buffer<int>(n_pl);

        // Get PL info
        int n = hts::bcf::get_info_int32(header, record, "DP", &buffer, &n_pl);
        BOOST_REQUIRE_GE(n,0);
        auto test = make_test_range(buffer.get(),buffer.get()+n);
        CHECK_EQUAL_RANGES(test, expected);
    }};

    test(make_data_url(vcfinfo), {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20});
}
