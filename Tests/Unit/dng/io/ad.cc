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
#include <dng/detail/varint.h>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace dng::io;
using namespace dng::pileup;
using namespace dng::utility;
using namespace std;

namespace dng { namespace io {
struct unittest_dng_io_ad {
    static int get_version_number(const dng::io::Ad &a) { return a.id_.version; }
    static std::string get_format_string(const dng::io::Ad &a) { return a.id_.name; }
};
}}
using dng::io::unittest_dng_io_ad;

// http://stackoverflow.com/a/673389
struct HexCharStruct {
    unsigned char c;
    HexCharStruct(unsigned char u) : c(u) { }
};

inline std::ostream& operator<<(std::ostream& o, const HexCharStruct& hs) {
    return (o << std::setw(2) << std::setfill('0') << std::hex << (int)hs.c << std::dec);
}

inline HexCharStruct hex(unsigned char u) {
    return HexCharStruct(u);
}

/*****************************************************************************
 Test the Varint format
 *****************************************************************************/

int ex_buf_size(uint64_t u) {
    int i = 1;
    for(; u >= 0x80; ++i, u >>= 7)
        /*noop*/;
    return i;
}

bool varint_convert(uint64_t u) {
    namespace varint = dng::detail::varint;
    stringbuf buffer;
    if(!varint::put(&buffer, u)) {
        cerr << "  Writing varint " << std::hex << std::showbase << u << " to a streambuf failed.\n";
        cerr << std::dec << std::noshowbase;
        return false;
    }
    int sz = ex_buf_size(u);
    if(buffer.str().length() != sz) {
        cerr << "  Writing varint " << std::hex << std::showbase << u << " to a streambuf failed.\n";
        cerr << "    Expected size: " << sz << " Actual size: " << buffer.str().length() << "\n";
        cerr << "    Buffer:";
        for(auto a : buffer.str()) {
            cerr << " " << hex(a);
        }
        cerr << "\n";
        cerr << std::dec << std::noshowbase;
        return false;
    }

    buffer.pubseekpos(0);
    auto aa = varint::get(&buffer);
    if(!aa.second) {
        cerr << "  Reading varint " << std::hex << std::showbase << u << " from a streambuf failed.\n";
        cerr << "    Buffer:";
        for(auto a : buffer.str()) {
            cerr << " " << hex(a);
        }
        cerr << "\n";
        cerr << std::dec << std::noshowbase;
        cerr << std::dec << std::noshowbase;
        return false;
    }
    if(aa.first != u) {
        cerr << "  Conversion of " << std::hex << std::showbase << u << " via streambuf failed.\n";
        cerr << "    Output: " << std::hex << std::showbase << aa.first << "\n";
        cerr << "    Buffer:";
        for(auto a : buffer.str()) {
            cerr << " " << hex(a);
        }
        cerr << "\n";
        cerr << std::dec << std::noshowbase;
        return false;
    }
    //cerr << "  Conversion of " << std::hex << std::showbase << u << " via streambuf success.\n";
    return true;
}

BOOST_AUTO_TEST_CASE(test_varint) {
    std::vector<uint64_t> numbers;
    {
        // use the lab's xorshift64 random generator
        uint64_t u = 15191868757011070976ULL;
        uint64_t w = 0x61C8864680B583EBLL;
        for(int i=0;i<10100;++i) {
            u ^= (u << 5); u ^= (u >> 15); u ^= (u << 27);
            w += 0x61C8864680B583EBLL;
            if(i < 100) {
                continue; //burnin 
            }
            numbers.push_back(u+(w^(w>>27)));
        }
    }

    for(int i=0;i<65536;++i) {
        BOOST_CHECK(varint_convert(i));
    }

    for(int i=0;i<64;++i) {
        BOOST_CHECK(varint_convert((1LL<<i)-1));
        BOOST_CHECK(varint_convert(1LL<<i));
        BOOST_CHECK(varint_convert((1LL<<i)+1));
    }

    for(int i=17;i<64;++i) {
        for(auto n : numbers) {

            BOOST_CHECK(varint_convert(n & ((1LL<<i)-1)));
        }
    }
}

/*****************************************************************************
 Test the Zig-Zag format
 *****************************************************************************/

int ex_bufz_size(int64_t n) {
    int i = 1;
    uint64_t u = (n < 0) ? -2*n-1 : 2*n;
    for(; u >= 0x80; ++i, u >>= 7)
        /*noop*/;
    return i;
}

bool zigzag_convert(int64_t n) {
    namespace varint = dng::detail::varint;
    stringbuf buffer;
    if(!varint::put_zig_zag(&buffer, n)) {
        cerr << "  Writing zig-zag varint " << std::hex << std::showbase << n << " to a streambuf failed.\n";
        cerr << std::dec << std::noshowbase;
        return false;
    }
    int sz = ex_bufz_size(n);
    if(buffer.str().length() != sz) {
        cerr << "  Writing zig-zag varint " << std::hex << std::showbase << n << " to a streambuf failed.\n";
        cerr << "    Expected size: " << sz << " Actual size: " << buffer.str().length() << "\n";
        cerr << "    Buffer:";
        for(auto a : buffer.str()) {
            cerr << " " << hex(a);
        }
        cerr << "\n";
        cerr << std::dec << std::noshowbase;
        return false;
    }

    buffer.pubseekpos(0);
    auto aa = varint::get_zig_zag(&buffer);
    if(!aa.second) {
        cerr << "  Reading zig-zag varint " << std::hex << std::showbase << n << " from a streambuf failed.\n";
        cerr << std::dec << std::noshowbase;
        return false;
    }
    if(aa.first != n) {
        cerr << "  Conversion of zig-zag " << std::hex << std::showbase << n << " via streambuf failed.\n";
        cerr << "    Output: " << std::hex << std::showbase << aa.first << "\n";
        cerr << "    Buffer:";
        for(auto a : buffer.str()) {
            cerr << " " << hex(a);
        }
        cerr << "\n";
        cerr << std::dec << std::noshowbase;
        return false;
    }
    return true;    
}

BOOST_AUTO_TEST_CASE(test_varint_zigzag) {
    std::vector<uint64_t> numbers;
    {
        // use the lab's xorshift64 random generator
        uint64_t u = 15191868757011070976ULL;
        uint64_t w = 0x61C8864680B583EBLL;
        for(int i=0;i<10100;++i) {
            u ^= (u << 5); u ^= (u >> 15); u ^= (u << 27);
            w += 0x61C8864680B583EBLL;
            if(i < 100) {
                continue; //burnin 
            }
            numbers.push_back(u+(w^(w>>27)));
        }
    }

    for(int i=0;i<65536;++i) {
        BOOST_CHECK(varint_convert(i));
        BOOST_CHECK(varint_convert(-i));
    }

    for(int i=0;i<64;++i) {
        int64_t j = (1LL << i);
        BOOST_CHECK(varint_convert(j-1));
        BOOST_CHECK(varint_convert(j));
        BOOST_CHECK(varint_convert(j+1));
        BOOST_CHECK(varint_convert(-j-1));
        BOOST_CHECK(varint_convert(-j));
        BOOST_CHECK(varint_convert(-j+1));
    }

    for(int i=17;i<=64;++i) {
        for(auto n : numbers) {
            BOOST_CHECK(varint_convert(n & ((1LL<<i)-1)));
            BOOST_CHECK(varint_convert(-(n & ((1LL<<i)-1))));
        }
    }
}


bool ad_write(std::vector<AlleleDepths> lines, std::vector<uint8_t> match) {
    std::stringstream buffer;
    Ad adfile("ad:", std::ios_base::out);
    adfile.AddContig("seq1",100);
    adfile.AddContig("seq2",1000, "");
    adfile.AddContig("seq3",10000, "M5:aaaaaaaa");
    adfile.AddLibrary("A","A");
    adfile.AddLibrary("B","B");
    adfile.Attach(buffer.rdbuf());
    for(auto line : lines) {
        if(adfile.Write(line) != 1) {
            std::cerr << "  Error: dng::io::Ad::Write failed\n";
            return false;
        }
    }    
    std::string result = buffer.str();

    if(result != std::string(match.begin(),match.end())) {
        std::cerr << "  Error: output does not match expectation\n";
        std::cerr << "  Result:   ";
        for(int i=0; i < result.size(); ++i) {
            std::cerr << hex(result[i]) << (i+1 < result.size() ? "," : "\n");
        }

        std::cerr << "  Expected: ";
        for(int i=0; i < match.size(); ++i) {
            std::cerr << hex(match[i]) << (i+1 < match.size() ? "," : "\n");
        }

        return false;
    }
    return true;
}

/*****************************************************************************
 Test the writing of a binary ad file
 *****************************************************************************/

BOOST_AUTO_TEST_CASE(test_ad_write) {
    BOOST_CHECK(ad_write({{make_location(0,0), 0, 2, {0,0}}}, {0x80,0x80,0x80,0x80,0x80,0x10,0x00,0x00}));
    BOOST_CHECK(ad_write({{make_location(1,9), 4, 2, {10,11,3,4}}}, {0x84,0x89,0x80,0x80,0x80,0x20,10,11,3,4}));
    BOOST_CHECK(ad_write({{make_location(2,100), 68, 2, {10,11,1000,4}}}, {0xC4,0xE4,0x80,0x80,0x80,0x30,10,11,0xE8,0x07,4}));

    BOOST_CHECK(ad_write({
        {make_location(1,0), 1, 2, {100,1001}},
        {make_location(1,1), 5, 2, {200,201,0,10}}, 
        {make_location(1,3), 5, 2, {200,201,0,10}}, 
        {make_location(2,3), 0, 2, {1,0}} 
    }, {
        0x81,0x80,0x80,0x80,0x80,0x20,    0x64, 0xE9,0x07,
                                 0x05,    0xC8,0x01, 0xBF,0x0C, 0,10,
                            0x85,0x01,    0x00, 0x00,0,10,
        0x80,0x83,0x80,0x80,0x80,0x30,    1, 0
    }));
}

/*****************************************************************************
 Test the writing and read of an ad file
 *****************************************************************************/
BOOST_AUTO_TEST_CASE(test_ad_write_and_read) {
    std::stringstream buffer;
    Ad adfile("ad:", std::ios_base::out | std::ios_base::in);
    adfile.Attach(buffer.rdbuf());

    adfile.AddContig("scaffold_1",1000);
    adfile.AddContig("scaffold_2",10000, "");
    adfile.AddContig("scaffold_3",100000, "M5:aaaaaaaa");
    adfile.AddLibrary("A","A");
    adfile.AddLibrary("B","B");

    std::vector<AlleleDepths> outdepths, indepths;
    int x=0;
    for(int i=0;i<128;++i) {
        AlleleDepths line;
        line.resize(i,2);
        line.location(0,10+(i+1));
        for(auto && d : line.data()) {
            d = x;
            x += 31;
        }
        outdepths.push_back(std::move(line));
    }
    for(int i=0;i<128;++i) {
        AlleleDepths line;
        line.resize(i,2);
        line.location(1,100+(i+1));
        for(auto && d : line.data()) {
            d = x;
            x += 31;
        }
        outdepths.push_back(std::move(line));
    }

    adfile.WriteHeader();
    for(auto && line : outdepths) {
        adfile.Write(line);
    }

    adfile.ReadHeader();

    BOOST_CHECK(unittest_dng_io_ad::get_version_number(adfile) == 0x0001);
    BOOST_CHECK(unittest_dng_io_ad::get_format_string(adfile) == "AD");

    BOOST_CHECK(adfile.contigs().size() == 3);
    BOOST_CHECK(adfile.contig(0).name == "scaffold_1");
    BOOST_CHECK(adfile.contig(1).name == "scaffold_2");
    BOOST_CHECK(adfile.contig(2).name == "scaffold_3");
    BOOST_CHECK(adfile.contig(0).length == 1000);
    BOOST_CHECK(adfile.contig(1).length == 10000);
    BOOST_CHECK(adfile.contig(2).length == 100000);

    BOOST_CHECK(adfile.libraries().size() == 2);
    BOOST_CHECK(adfile.library(0).name == "A");
    BOOST_CHECK(adfile.library(1).name == "B");
    BOOST_CHECK(adfile.library(0).sample == "A");
    BOOST_CHECK(adfile.library(1).sample == "B");

    AlleleDepths input;
    while(adfile.Read(&input)) {
        indepths.push_back(input);
    }
    BOOST_CHECK(outdepths.size() == indepths.size());
    for(int i=0;i<128;++i) {
        const AlleleDepths &a = outdepths[i];
        const AlleleDepths &b = indepths[i];
        BOOST_CHECK(a.location() == b.location());
        BOOST_CHECK(a.color() == b.color());
        BOOST_CHECK(a.data() == b.data());        
    }
}

/*****************************************************************************
 Test the writing of a tad file
 *****************************************************************************/

bool tad_write(std::vector<AlleleDepths> lines, std::string text) {
    std::stringstream buffer;
    Ad adfile("tad:", std::ios_base::out);
    adfile.Attach(buffer.rdbuf());

    adfile.AddContig("seq1",100);
    adfile.AddContig("seq2",1000, "");
    adfile.AddContig("seq3",10000, "M5:aaaaaaaa");
    adfile.AddLibrary("A","A");
    adfile.AddLibrary("B","B");

    for(auto line : lines) {
        if(adfile.Write(line) != 1) {
            std::cerr << "  Error: Ad::Write failed\n";
            return false;
        }
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

BOOST_AUTO_TEST_CASE(test_tad_write) {
    BOOST_CHECK(tad_write({{make_location(0,0), 0, 2, {0,0}}}, "seq1\t1\tA\t0\t0\n"));
    BOOST_CHECK(tad_write({{make_location(1,9), 4, 2, {10,11,3,4}}}, "seq2\t10\tAC\t10,3\t11,4\n"));
    BOOST_CHECK(tad_write({{make_location(2,100), 68, 2, {10,11,3,4}}}, "seq3\t101\tNAC\t10,3\t11,4\n"));

    BOOST_CHECK(tad_write({
        {make_location(0,0), 1, 2, {100,1001}}, {make_location(0,1), 5, 2, {200,201,0,10}} 
    }, "seq1\t1\tC\t100\t1001\n" "seq1\t2\tAG\t200,0\t201,10\n"));
}

/*****************************************************************************
 Test the reading of a tad file
 *****************************************************************************/

const char tad_test1[] = 
    "@ID\tFF:TAD\tVN:0.1\n"
    "@SQ\tSN:scaffold_1\tLN:100\n"
    "@SQ\tSN:scaffold_2\tLN:200\tM5:aaaaaaaa\n"
    "@SQ\tSN:scaffold_3\tLN:300\tM5:aaaaaaab\tUR:blah\n"
    "@AD\tID:A\tSM:A\n"
    "@AD\tID:B\tSM:B\tLB:B\tRG:B\n"
    "@CO\tThis is a comment that is ignored\n"
    "\n"
    "scaffold_1\t1\tA\t10\t9\n"
    "scaffold_1\t2\tC\t8\t7\n"
    "scaffold_1\t3\tGC\t6,0\t0,2\n"
    "scaffold_2\t1\tT\t4\t0\n"
    "scaffold_3\t1\tntgca\t0,0,0,1\t4,0,0,0\n"
;

BOOST_AUTO_TEST_CASE(test_tad_read) {
    std::stringstream buffer(tad_test1);
    Ad adfile("tad:", std::ios_base::in);
    adfile.Attach(buffer.rdbuf());
    adfile.ReadHeader();

    typedef std::vector<int> V;
    typedef std::vector<std::string> S;


    BOOST_CHECK(unittest_dng_io_ad::get_version_number(adfile) == 0x0001);
    BOOST_CHECK(unittest_dng_io_ad::get_format_string(adfile) == "TAD");

    BOOST_CHECK(adfile.contigs().size() == 3);
    BOOST_CHECK(adfile.contig(0).name == "scaffold_1");
    BOOST_CHECK(adfile.contig(1).name == "scaffold_2");
    BOOST_CHECK(adfile.contig(2).name == "scaffold_3");
    BOOST_CHECK(adfile.contig(0).length == 100);
    BOOST_CHECK(adfile.contig(1).length == 200);
    BOOST_CHECK(adfile.contig(2).length == 300);
    BOOST_CHECK(adfile.contig(0).attributes.empty());
    BOOST_CHECK(adfile.contig(1).attributes == S({"M5:aaaaaaaa"}));
    BOOST_CHECK(adfile.contig(2).attributes == S({"M5:aaaaaaab","UR:blah"}));

    AlleleDepths depths;
    adfile.Read(&depths);
    BOOST_CHECK(depths.location() == make_location(0,0));
    BOOST_CHECK(depths.color() == 0);
    BOOST_CHECK(depths.data() == V({10,9}));

    adfile.Read(&depths);
    BOOST_CHECK(depths.location() == make_location(0,1));
    BOOST_CHECK(depths.color() == 1);
    BOOST_CHECK(depths.data() == V({8,7}));

    adfile.Read(&depths);
    BOOST_CHECK(depths.location() == make_location(0,2));
    BOOST_CHECK(depths.color() == 11);
    BOOST_CHECK(depths.data() == V({6,0,0,2}));

    adfile.Read(&depths);
    BOOST_CHECK(depths.location() == make_location(1,0));
    BOOST_CHECK(depths.color() == 3);
    BOOST_CHECK(depths.data() == V({4,0}));    

    adfile.Read(&depths);
    BOOST_CHECK(depths.location() == make_location(2,0));
    BOOST_CHECK(depths.color() == 127);
    BOOST_CHECK(depths.data() == V({0,4,0,0,0,0,1,0}));
}

bool string_equal(const std::string &a, const std::string &b) {
    if(a == b) {
        return true;
    }
    std::cerr << "  Error: output does not match expectation\n";
    std::cerr << "  Result:\n" << a << "\n";
    std::cerr << "  Expected:\n" << b << "\n";
    return false;
}

const char tad_test2[] = 
    "@ID\tFF:TAD\tVN:0.1\n"
    "@SQ\tSN:scaffold_1\tLN:100\n"
    "@SQ\tSN:scaffold_2\tLN:200\tM5:aaaaaaaa\n"
    "@SQ\tSN:scaffold_3\tLN:300\tM5:aaaaaaab\tUR:blah\n"
    "@AD\tID:A\tSM:A\n"
    "@AD\tID:B\tSM:B\tLB:B\tRG:B\n"
    "@CO\tThis is a comment that is ignored\n"
    "scaffold_1\t1\tA\t10\t9\n"
    "scaffold_1\t2\tC\t8\t7\n"
    "scaffold_1\t3\tGC\t6,0\t0,2\n"
    "scaffold_2\t1\tT\t4\t0\n"
    "scaffold_3\t1\tNTGCA\t0,0,0,1\t4,0,0,0\n"
;

BOOST_AUTO_TEST_CASE(test_tad_read_and_write) {
    std::stringstream buffer(tad_test2);
    std::stringstream buffer_out;
    Ad adfile("tad:", std::ios_base::in | std::ios_base::out);
    adfile.Attach(buffer.rdbuf());
    adfile.ReadHeader();
    std::vector<AlleleDepths> depths_store;
    AlleleDepths depths;
    while(adfile.Read(&depths)) {
        depths_store.push_back(depths);
    }
    adfile.Attach(buffer_out.rdbuf());
    adfile.WriteHeader();
    for(auto &&a : depths_store) {
        adfile.Write(a);
    }
    BOOST_CHECK(string_equal(buffer_out.str(), tad_test2));
}


