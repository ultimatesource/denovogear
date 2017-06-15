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

#include <fstream>

using namespace dng::regions;

namespace dng { namespace regions {
bool operator==(const parsed_range_t &a, const parsed_range_t &b) {
    return std::tie(a.target,a.beg,a.end) == std::tie(b.target,b.beg,b.end);
}
}}

namespace hts { namespace bam {
bool operator==(const region_t &a, const region_t &b) {
    return std::tie(a.tid,a.beg,a.end) == std::tie(b.tid,b.beg,b.end);
}
}}

bool region_parsing(std::string input, parsed_ranges_t b ) {
    std::cerr << "Parsing region '" << input << "': ";
    auto result = parse_ranges(input);
    if(result.second == false) {
        std::cerr << "FAIL\n";
        return false;
    } else {
        std::cerr << "SUCCESS\n";
    }
    if(result.first == b) {
        return true;
    }
	std::cerr << "  Error: output does not match expectation\n";
    for(auto &&aa : result.first) {
        std::cerr << "  result: " << aa.target << "\t" << aa.beg << "\t" << aa.end << "\n";
    }
    return false;
}

bool region_parsing_expect_fail(std::string input) {
    std::cerr << "Parsing region '" << input << "': ";
    auto result = parse_ranges(input);
    if(result.second == false) {
        std::cerr << "FAIL (expected)\n";
        return true;
    } else {
        std::cerr << "SUCCESS (unexpected)\n";
    }
    for(auto &&aa : result.first) {
        std::cerr << "  result: " << aa.target << "\t" << aa.beg << "\t" << aa.end << "\n";
    }
    return false;
}

bool region_bam_parsing(std::string input, hts::bam::File &file, hts::bam::regions_t b ) {
    hts::bam::regions_t a;
    try {
        a = bam_parse_region(input, file);
    } catch(std::exception &e) {
        std::cerr << e.what() << "\n";
        return false;
    }
    if(a == b) {
        return true;
    }
    for(auto &&aa : a) {
        std::cerr << "    Parsing result: " << aa.tid << "\t" << aa.beg << "\t" << aa.end << "\n";
    }

    return false;
}

bool region_bam_parsing_expect_fail(std::string input, hts::bam::File &file) {
    try {
        bam_parse_region(input, file);
    } catch(std::exception &e) {
        std::cerr << "    Expected exception: " << e.what() << "\n";
        return true;
    }
    return false;
}


const char bamtext[]=
"@HD\tVN:1.4\tGO:none\tSO:coordinate\n"
"@SQ\tSN:1\tLN:249250621\n"
"@SQ\tSN:2\tLN:243199373\n"
;

BOOST_AUTO_TEST_CASE(test_region_parsing) {
    BOOST_CHECK(region_parsing_expect_fail("chr1:-"));
    BOOST_CHECK(region_parsing_expect_fail("chr1:"));
    BOOST_CHECK(region_parsing_expect_fail("chr1:1-1a1"));
    
    BOOST_CHECK(region_parsing("chr1:10-1000", {{"chr1",10,1000}}));
    BOOST_CHECK(region_parsing("chr1:10-1e3", {{"chr1",10,1000}}));
    BOOST_CHECK(region_parsing("chr1:10-1.2e3", {{"chr1",10,1200}}));
    BOOST_CHECK(region_parsing("chr1 : 10 - 1000", {{"chr1",10,1000}}));
    BOOST_CHECK(region_parsing("chr1:10-",     {{"chr1",10,INT_MAX}}));
    BOOST_CHECK(region_parsing("chr1:-1000",   {{"chr1",1,1000}}));
    BOOST_CHECK(region_parsing("chr1:1000",    {{"chr1",1000,1000}}));
    BOOST_CHECK(region_parsing("chr1:1,000",    {{"chr1",1000,1000}}));
    BOOST_CHECK(region_parsing("chr1:1.0e3",    {{"chr1",1000,1000}}));
    BOOST_CHECK(region_parsing("chr1",         {{"chr1",1,INT_MAX}}));
    BOOST_CHECK(region_parsing("chr1:100,001 -1,001,001",    {{"chr1",100001,1001001}}));

    BOOST_CHECK(region_parsing("chr1:10+1000", {{"chr1",10,1009}}));

    BOOST_CHECK(region_parsing("chr1:10-1000 3:60-100", {{"chr1",10,1000},{"3",60,100}}));
    BOOST_CHECK(region_parsing("chr1:10-1000  3", {{"chr1",10,1000},{"3",1,INT_MAX}}));
    BOOST_CHECK(region_parsing("\nchr1\t:\n10\f-\v1000\t\n \f\v3 ", {{"chr1",10,1000},{"3",1,INT_MAX}}));
    
    // Open Read our test data from Memory
    std::string bamstr = "data:";
    if(hts::version() >= 10400 ) {
        // data: url format changed in htslib 1.4
        bamstr += ",";
    }
    bamstr += bamtext;

    hts::bam::File bamfile(bamstr.c_str(), "r");
    // Sanity checks.  If these fail, then other tests will as well
    BOOST_CHECK(bamfile.is_open() && bamfile.header() != nullptr);
    BOOST_CHECK(bamfile.TargetNameToID("xyz") == -1);
    BOOST_CHECK(bamfile.TargetNameToID("1") == 0);
    BOOST_CHECK(bamfile.TargetNameToID("2") == 1);
    // Check normal 
    BOOST_CHECK(region_bam_parsing("1:1-1000", bamfile, {{0,0,1000}}));
    BOOST_CHECK(region_bam_parsing("1:1-1000 2:1", bamfile, {{0,0,1000},{1,0,1}}));
    // Check sorting
    BOOST_CHECK(region_bam_parsing("2:100-120 2:21-30 2:1-20 1:1-5 1:7-10", bamfile, {{0,0,5},{0,6,10},{1,0,30},{1,99,120}}));
    // Check merging
    BOOST_CHECK(region_bam_parsing("1:1-1000 1:900-2000", bamfile, {{0,0,2000}}));
    BOOST_CHECK(region_bam_parsing("2:10-200 2:100-300 1:1-1000 2:200-400", bamfile, {{0,0,1000},{1,9,400}}));
    // Check expected failures
    BOOST_CHECK(region_bam_parsing_expect_fail("xyz", bamfile));
    BOOST_CHECK(region_bam_parsing_expect_fail("1:1000-900", bamfile));
    BOOST_CHECK(region_bam_parsing_expect_fail("1:1000-999", bamfile));
    // Check empty file
    BOOST_CHECK(region_bam_parsing_expect_fail("", bamfile));
}

bool bed_parsing(std::string input, hts::bam::File &file, hts::bam::regions_t b ) {
    hts::bam::regions_t a;

    dng::detail::AutoTempFile output;

    if(!output.file.is_open()) {
        std::cerr << "    Failed to open temporary file '" + output.path.string() + "'\n";
        return false;
    }
    output.file.write(input.c_str(),input.size());
    output.file.flush();

    try {
        a = bam_parse_bed(output.path.string(), file);
    } catch(std::exception &e) {
        std::cerr << e.what() << "\n";
        return false;
    }
    if(a == b) {
        return true;
    }
    for(auto &&aa : a) {
        std::cerr << "    Parsing result: " << aa.tid << "\t" << aa.beg << "\t" << aa.end << "\n";
    }

    return false;
}

BOOST_AUTO_TEST_CASE(test_bed_parsing) {
    // Open Read our test data from Memory
    std::string bamstr = "data:";
    if(hts::version() >= 10400 ) {
        // data: url format changed in htslib 1.4
        bamstr += ",";
    }
    bamstr += bamtext;

    hts::bam::File bamfile(bamstr.c_str(), "r");
    // Sanity checks.  If these fail, then other tests will as well
    BOOST_CHECK(bamfile.is_open() && bamfile.header() != nullptr);
    BOOST_CHECK(bamfile.TargetNameToID("xyz") == -1);
    BOOST_CHECK(bamfile.TargetNameToID("1") == 0);
    BOOST_CHECK(bamfile.TargetNameToID("2") == 1);
    // Check normal 
    BOOST_CHECK(bed_parsing("1\t0\t1000\n", bamfile, {{0,0,1000}}));
    BOOST_CHECK(bed_parsing("1\t0\t1000\n2\t0\t1", bamfile, {{0,0,1000},{1,0,1}}));
    // Check sorting
    BOOST_CHECK(bed_parsing("2\t99\t120\n2\t20\t30\n2\t0\t20\n1\t0\t5\n1\t6\t10", bamfile, {{0,0,5},{0,6,10},{1,0,30},{1,99,120}}));
    // Check merging
    BOOST_CHECK(bed_parsing("1\t0\t1000\n1\t899\t2000", bamfile, {{0,0,2000}}));
    BOOST_CHECK(bed_parsing("2\t9\t200\n2\t99\t300\n1\t0\t1000\n2\t199\t400", bamfile, {{0,0,1000},{1,9,400}}));
    // Check Comments and blank lines
    BOOST_CHECK(bed_parsing("#comment\n\n\n#comment\ttest\n#comment\ttest\tanothertest\n1\t0\t1000\n", bamfile, {{0,0,1000}}));
    // Check empty file
    BOOST_CHECK(bed_parsing("", bamfile, {}));

}
