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

BOOST_AUTO_TEST_CASE(test_query_length) {
    auto test = [](cigar_t cigar, std::size_t expected_result) -> void {
    BOOST_TEST_CONTEXT("cigar= (" << *cigar.first << "," << *cigar.second << ")") {
	auto result = query_length(cigar);
	BOOST_CHECK_EQUAL(result, expected_result);
    }};
    
    std::vector<uint32_t> v = {4000, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {81, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {82, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
    v = {83, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
    v = {84, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {85, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
    v = {86, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
    v = {87, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {88, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {89, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
    v = {0, 0};
    test(std::make_pair(&v.front(),&v.back()),0);
}

using dng::cigar::target_length;

BOOST_AUTO_TEST_CASE(test_target_length) {
    auto test = [](cigar_t cigar, std::size_t expected_result) -> void {
    BOOST_TEST_CONTEXT("cigar= (" << *cigar.first << "," << *cigar.second << ")") {
        auto result = target_length(cigar);
        BOOST_CHECK_EQUAL(result, expected_result);
    }};

    std::vector<uint32_t> v = {4000, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {81, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
    v = {82, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {83, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {84, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
    v = {85, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
    v = {86, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
    v = {87, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {88, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),250);
    v = {89, 3920, 0000};
    test(std::make_pair(&v.front(),&v.back()),245);
}

using dng::cigar::query_pos;

BOOST_AUTO_TEST_CASE(test_query_pos) {
    auto test = [](uint64_t q, uint64_t expected_result) -> void {
    BOOST_TEST_CONTEXT("query= " << q) {
	auto result = query_pos(q);
	BOOST_CHECK_EQUAL(result, expected_result);
    }};
    
    test(126385667, 63192833);
    test(126385666, 63192833);
    test(126385420, 63192710);
}

using dng::cigar::query_del;

BOOST_AUTO_TEST_CASE(test_query_del) {
    auto test = [](uint64_t q, uint64_t expected_result) -> void {
    BOOST_TEST_CONTEXT("query= " << q) {
        auto result = query_del(q);
        BOOST_CHECK_EQUAL(result, expected_result);
    }};

    test(126385667, 1);
    test(126385666, 0);
    test(126385420, 0);
}

using dng::cigar::target_to_query;

BOOST_AUTO_TEST_CASE(test_target_to_query) {
    auto test = [](uint64_t target, uint64_t beg, cigar_t cigar, uint64_t expected_result) -> void {
    BOOST_TEST_CONTEXT("target= " << target << ", beg=" << beg) {
        auto result = target_to_query(target, beg, cigar);
        BOOST_CHECK_EQUAL(result, expected_result);
    }};

    std::vector<uint32_t> v = {4000,0000};
    test(126385670, 126385665, std::make_pair(&v.front(), &v.back()), 10);
    test(126385670, 126385670, std::make_pair(&v.front(), &v.back()), 0);
    test(126385650, 60192825, std::make_pair(&v.front(), &v.back()), 499);
    test(126385680, 126385685, std::make_pair(&v.front(), &v.back()), -1); 
    v = {81, 3920, 0000};
    test(126385670, 126385665, std::make_pair(&v.front(), &v.back()), 20);
    test(126385670, 126385670, std::make_pair(&v.front(), &v.back()), 10);
    test(126385650, 60192825, std::make_pair(&v.front(), &v.back()), 499);
    test(126385680, 126385685, std::make_pair(&v.front(), &v.back()), -1);
}
