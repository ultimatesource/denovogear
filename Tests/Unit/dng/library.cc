/*
 * Copyright (c) 2017 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *	     Juan J. Garcia Mesa <jgarc111@asu.edu>
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

#define BOOST_TEST_MODULE dng::library

#include <dng/library.h>

#include "../testing.h"

using namespace dng;

using dng::trim_label_prefix;

BOOST_AUTO_TEST_CASE(test_remove_prefixes_labels) {
    auto test = [](const std::string label, const std::string  expected_result) -> void {
    BOOST_TEST_CONTEXT("label= " << label) {
	auto result = trim_label_prefix(label);
	BOOST_CHECK_EQUAL(result, expected_result);
    }};

    test("LB/NA12891", "NA12891");
    test("LB/NA12878/Solexa-135852", "NA12878/Solexa-135852");
    test("LB/NA12891/H03N7", "NA12891/H03N7");
    test("LB//NA12891", "/NA12891");
    test("GL/NA12891", "NA12891");
    test("GL/1", "1");
    test("SM/NA12891", "NA12891");
    test("NA12892", "NA12892");
    test("/NA12892", "/NA12892");
    test(" LB/NA128921", " LB/NA128921");
    test(" ", " ");
}

using dng::trim_label_prefix_only_libraries;

BOOST_AUTO_TEST_CASE(test_prefix_ony_libraries) {
    auto test = [](const std::string label, const std::string expected_result) -> void {
    BOOST_TEST_CONTEXT("label= " << label) {
	auto result = trim_label_prefix_only_libraries(label);
	BOOST_CHECK_EQUAL(result, expected_result);
    }};

    test("LB/NA12891", "NA12891");
    test("LB/NA12878/Solexa-135852", "NA12878/Solexa-135852");
    test("LB/NA12891/H03N7", "NA12891/H03N7");    
    test("LB//NA12891", "/NA12891");
    test(" LB/NA12891", " LB/NA12891");
    test("GL/NA12891", "");
    test("GL/1", "");
    test("SM/NA12891", "");
    test("NA12892", "NA12892");
    test("/NA12892", "/NA12892");
    test(" LB/NA12891", " LB/NA12891");
    test(" ", " ");
}
