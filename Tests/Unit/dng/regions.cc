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
#include <fstream>

using namespace dng::regions;

namespace dng { namespace regions { namespace detail {
bool operator==(const raw_parsed_region_t &a, const raw_parsed_region_t &b) {
    return std::tie(a.target,a.beg,a.end) == std::tie(b.target,b.beg,b.end);
}
}}}

namespace hts { namespace bam {
bool operator==(const region_t &a, const region_t &b) {
    return std::tie(a.tid,a.beg,a.end) == std::tie(b.tid,b.beg,b.end);
}
}}

bool region_parsing(std::string input, detail::raw_parsed_regions_t b ) {
    std::cerr << "Parsing region '" << input << "': ";
    auto result = detail::parse_regions(input);
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
    auto result = detail::parse_regions(input);
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

const char bamtext[]= "data:"
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
    hts::bam::File bamfile(bamtext, "r");
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
}

bool region_bed_parsing(std::string input, hts::bam::File &file, hts::bam::regions_t b ) {
    return false; // figure out how to autodelete temp file

    hts::bam::regions_t a;

    std::string path = std::tmpnam(nullptr);
    std::ofstream output{path};

    if(!output) {
        std::cerr << "    Failed to open temporary file.\n";
        return false;
    }
    output.write(input.c_str(),input.size());

    try {
        a = bam_parse_bed(path, file);
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
