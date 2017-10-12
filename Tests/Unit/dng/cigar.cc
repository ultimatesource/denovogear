/*
 * Copyright (c) 2017 Reed A. Cartwright
 * Authors:  Juan J. Garcia Mesa <jgarc111@asu.edu>
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

#define BOOST_TEST_MODULE dng::cigar

#include <vector>

#include <dng/cigar.h>

#include "../testing.h"

using namespace dng;

using dng::cigar::query_length;
using hts::bam::cigar_t;

std::vector<uint32_t> cigar_v(int first, int second){
    std::vector<uint32_t> v;
    while(first<second) {
	v.push_back(first);
	++first;
    }
    return v;
}

BOOST_AUTO_TEST_CASE(test_query_length) {
    auto test = [](cigar_t cigar, std::size_t expected_result) -> void {
    BOOST_TEST_CONTEXT("cigar= (" << *cigar.first << "," << *cigar.second << ")") {
	auto result = query_length(cigar);
	BOOST_CHECK_EQUAL(result, expected_result);
    }};
    
    std::vector<uint32_t> v = cigar_v(39416484, 39416488);
    cigar_t c = std::make_pair(&v.front(), &v.back());
    test(c, 2463530);
}

using dng::cigar::target_to_query;

BOOST_AUTO_TEST_CASE(test_target_to_query) {
    auto test = [](uint64_t target, uint64_t beg, cigar_t cigar, uint64_t expected_result) -> void {
    BOOST_TEST_CONTEXT("cigar= (" << *cigar.first << "," << *cigar.second << ")") {
	auto result = target_to_query(target, beg, cigar);
	BOOST_CHECK_EQUAL(result, expected_result);
    }};

    std::vector<uint32_t> v = cigar_v(39416484, 39416489);
    cigar_t c = std::make_pair(&v.front(), &v.back());
    test(0, 1, c, -1);
    test(2, 1, c, 4927062);
    test(100000000, 1, c, 9854119);
}

using dng::cigar::target_length;

BOOST_AUTO_TEST_CASE(test_target_length) {
    auto test = [](cigar_t cigar, std::size_t expected_result) -> void {
    BOOST_TEST_CONTEXT("cigar= (" << *cigar.first << "," << *cigar.second << ")") {
	auto result = target_length(cigar);
	BOOST_CHECK_EQUAL(result, expected_result);
    }};
    
    std::vector<uint32_t> v = cigar_v(39416484, 39416489);
    cigar_t c = std::make_pair(&v.front(), &v.back()); 
    test(c, 2463530);
        
}

using dng::cigar::query_pos;

BOOST_AUTO_TEST_CASE(test_query_pos) {
    auto test = [](uint64_t q, uint64_t expected_result) -> void {
    BOOST_TEST_CONTEXT("query= " << q) {
	auto result = query_pos(q);
	BOOST_CHECK_EQUAL(result, expected_result);
    }};
    
    test(2, 1);
    test(0, 0);
}

using dng::cigar::query_del;

BOOST_AUTO_TEST_CASE(test_query_del) {
    auto test = [](uint64_t q, uint64_t expected_result) -> void {
    BOOST_TEST_CONTEXT("query= " << q) {
        auto result = query_del(q);
        BOOST_CHECK_EQUAL(result, expected_result);
    }};

    test(2, 0);
    test(1, 1);
}
