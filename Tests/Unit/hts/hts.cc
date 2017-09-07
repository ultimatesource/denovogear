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

#define BOOST_TEST_MODULE hts::hts

#include <dng/hts/hts.h>
#include "../testing.h"

using namespace hts;
using hts::version_parse;

BOOST_AUTO_TEST_CASE(test_hts_version) {
    auto test = [](const char* version_char, unsigned long version_int) ->  void {    
    BOOST_TEST_CONTEXT("version_string='" << std::string(version_char) << "'"){
	auto result = version_parse(version_char);
	BOOST_CHECK_EQUAL(result, version_int);
    }};

    test("1.5-16-g7cdd574",10516);
    test("1.4.1",10401);
    test("1.4",10400);
    test("001.4.01",10401);
    test("1.4-1",10401);
    test("0",0);
    test("htslib version 1.5",0);
    test("1,5",0);
    test("1 . 4",0);
    test(" ",0);
}
