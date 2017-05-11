/*
 * Copyright (c) 2017 Reed A. Cartwright
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
#define BOOST_TEST_MODULE dng::io::bam

#include <dng/detail/unit_test.h>

#include <dng/io/bam.h>

#include <vector>

namespace dng { namespace io {
struct unittest_dng_io_bam {
    GETTERF(BamPileup, ParseHeader);
};
}}

using namespace dng;
using namespace dng::io;
using namespace std;

const char header1[] =
    "@RG\tID:Mom\tLB:Mom\tSM:Mom\n"
    "@RG\tID:Dad1\tLB:Dad\tSM:Dad\n"
    "@RG\tID:Dad2\tLB:Dad\tSM:Dad\n"
    "@RG\tID:Eve1\tLB:Eve1\tSM:Eve\n"
    "@RG\tID:Eve2\tLB:Eve2\tSM:Eve\n"
;

BOOST_AUTO_TEST_CASE(test_header_parser) {
    {
        BamPileup pileup;
        vector<string> expected_names = {"Mom","Dad", "Eve1", "Eve2"};
        vector<string> expected_samples = {"Mom","Dad", "Eve", "Eve"};

        unittest_dng_io_bam::ParseHeader(pileup, header1);
        BOOST_CHECK_EQUAL(pileup.num_libraries(), 4);
        CHECK_EQUAL_RANGES(pileup.libraries().names, expected_names);
        CHECK_EQUAL_RANGES(pileup.libraries().samples, expected_samples);
    }
}
