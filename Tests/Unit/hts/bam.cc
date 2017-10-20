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

#define BOOST_TEST_MODULE hts::bam

#include <dng/hts/bam.h>

#include "../testing.h"
#include <dng/hts/hts.h>

using namespace hts;
using namespace hts::bam;
using hts::detail::make_data_url;

const char samwithflags[] = 
    "@HD\tVN:1.5\tSO:coordinate\n"
    "@SQ\tSN:ref\tLN:10\n"
    "r001\t66\tref\t1\t40\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n"
    "r002\t4\tref\t1\t40\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n"  // BAM_FUNMAP removed
    "r003\t258\tref\t1\t40\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n" // BAM_FSECONDARY removed
    "r004\t514\tref\t1\t40\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n" // BAM_FQCFAIL removed
    "r005\t1026\tref\t1\t40\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n" // BAM_FDUP removed
    "r006\t2050\tref\t1\t40\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n" // BAM_FSUPPLEMENTARY
    "r007\t82\tref\t1\t40\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n"    
;

BOOST_AUTO_TEST_CASE(test_file_read) {
    std::string bamstr = make_data_url(samwithflags);
    hts::bam::File bamfile(bamstr.c_str(), "r");
    Alignment aln;
    BOOST_REQUIRE_GE(bamfile.Read(&aln),0);
    BOOST_CHECK_EQUAL(aln.qname(), "r001");
    BOOST_REQUIRE_GE(bamfile.Read(&aln),0);
    BOOST_CHECK_EQUAL(aln.qname(), "r007");
    BOOST_CHECK_LE(bamfile.Read(&aln),0);
}
