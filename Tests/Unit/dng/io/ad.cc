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
#define BOOST_TEST_MODULE dng::io::ad

#include <dng/io/ad.h>
#include <sstream>

using namespace dng::io;
using namespace dng::pileup;
using namespace dng::utility;

bool tad_write(AlleleDepths line, std::string text) {
    std::stringstream buffer;
    Ad adfile("tad:", std::ios_base::out);
    adfile.contigs({"seq1", "seq2", "seq3"});
    adfile.Attach(buffer.rdbuf());
    if(adfile.Write(line) != 0) {
        std::cerr << "  Error: adfile.Write failed\n";
        return false;
    }
    std::string ss = buffer.str();
    if(ss != text) {
        std::cerr << "  Error: output does not match expectation\n";
        std::cerr << "  Result:   " << ss << "\n";
        std::cerr << "  Expected: " << text << "\n";
        return false;
    }
    return true;
}

BOOST_AUTO_TEST_CASE(test_ad_write) {
    BOOST_CHECK(tad_write({make_location(0,0), 0, 2, {0,0}}, "seq1\t1\tA\t0\t0\n"));
    BOOST_CHECK(tad_write({make_location(1,9), 4, 2, {10,11,3,4}}, "seq2\t10\tAC\t10,3\t11,4\n"));
    BOOST_CHECK(tad_write({make_location(2,100), 68, 2, {10,11,3,4}}, "seq3\t101\tNAC\t10,3\t11,4\n"));
}

