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
#include <iostream>
#include <iomanip>

using namespace dng::io;
using namespace dng::pileup;
using namespace dng::utility;

// http://stackoverflow.com/a/673389
struct HexCharStruct
{
  unsigned char c;
  HexCharStruct(unsigned char u) : c(u) { }
};

inline std::ostream& operator<<(std::ostream& o, const HexCharStruct& hs)
{
  return (o << std::setw(2) << std::setfill('0') << std::hex << (int)hs.c);
}

inline HexCharStruct hex(unsigned char u)
{
  return HexCharStruct(u);
}

/*****************************************************************************
 Test the NTF8 format
 *****************************************************************************/

int ntf8_put32(uint32_t n, char *out, size_t count);
int ntf8_put64(uint64_t n, char *out, size_t count);
int ntf8_get32(const char *in, size_t count, uint32_t *r);
int ntf8_get64(const char *in, size_t count, uint64_t *r);

bool ntf8_convert32(uint32_t n) {
    char buffer[5];
    int sz1 = ntf8_put32(n,buffer,5);
    uint32_t u;
    int sz2 = ntf8_get32(buffer,5,&u);
    if(sz1 != sz2 || u != n) {
        std::cerr << "  NTF8 32 conversion of " << n << "failed.\n";
        std::cerr << "    Output: " << u << ", size written = " << sz1 << ", size read = " << sz2 << ".\n";
        return false;
    }

    //std::cerr << "  NTF8 32 conversion of " << n << ", size read/written = " << sz1 << " successful.\n";
    return true;
}

bool ntf8_convert64(uint64_t n) {
    char buffer[9];
    int sz1 = ntf8_put64(n,buffer,9);
    uint64_t u;
    int sz2 = ntf8_get64(buffer,9,&u);
    if(sz1 != sz2 || u != n) {
        std::cerr << "  NTF8 64 conversion of " << n << "failed.\n";
        std::cerr << "    Output: " << u << ", size written = " << sz1 << ", size read = " << sz2 << ".\n";
        return false;
    }

    //std::cerr << "  NTF8 64 conversion of " << n << ", size read/written = " << sz1 << " successful.\n";
    return true;    
}


BOOST_AUTO_TEST_CASE(test_ntf8) {
    std::vector<uint64_t> numbers;
    {
        // user the lab's xorshift64 random generator
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
        BOOST_CHECK(ntf8_convert32(i));
    }

    for(int i=0;i<32;++i) {
        BOOST_CHECK(ntf8_convert32((1<<i)-1));
        BOOST_CHECK(ntf8_convert32(1<<i));
        BOOST_CHECK(ntf8_convert32((1<<i)+1));
    }

    for(int i=17;i<=32;++i) {
        for(auto n : numbers) {
            BOOST_CHECK(ntf8_convert32(n & ((1LL<<i)-1)));
        }
    }

    for(int i=0;i<65536;++i) {
        BOOST_CHECK(ntf8_convert64(i));
    }

    for(int i=0;i<64;++i) {
        BOOST_CHECK(ntf8_convert64((1LL<<i)-1));
        BOOST_CHECK(ntf8_convert64(1LL<<i));
        BOOST_CHECK(ntf8_convert64((1LL<<i)+1));
    }

    for(int i=17;i<=64;++i) {
        for(auto n : numbers) {
            BOOST_CHECK(ntf8_convert32(n & ((1LL<<i)-1)));
        }
    }

}

bool ad_write(std::vector<AlleleDepths> lines, std::vector<uint8_t> match) {
    std::stringstream buffer;
    Ad adfile("ad:", std::ios_base::out);
    adfile.AddContigs("seq1",100);
    adfile.AddContigs("seq2",1000, "");
    adfile.AddContigs("seq3",10000, "M5:aaaaaaaa");
    adfile.Attach(buffer.rdbuf());
    for(auto line : lines) {
        if(adfile.Write(line) != 0) {
            std::cerr << "  Error: adfile.Write failed\n";
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
    BOOST_CHECK(ad_write({{make_location(0,0), 0, 2, {0,0}}}, {0xF8,0x80,0x00,0x00,0x00,0x00,0x00,0x00}));
    BOOST_CHECK(ad_write({{make_location(1,9), 4, 2, {10,11,3,4}}}, {0xF9,0x00,0x00,0x00,0x04,0x84,10,11,3,4}));
    BOOST_CHECK(ad_write({{make_location(2,100), 68, 2, {10,11,1000,4}}}, {0xF9,0x80,0x00,0x00,0x32,0x44,10,11,0x83,0xE8,4}));

    BOOST_CHECK(ad_write({
        {make_location(1,0), 1, 2, {100,1001}},
        {make_location(1,1), 5, 2, {200,201,0,10}}, 
        {make_location(1,3), 5, 2, {200,201,0,10}}, 
        {make_location(2,3), 0, 2, {1,0}} 
    }, {
        0xF9,0x00,0x00,0x00,0x00,0x01,    0x64, 0x83,0xE9,
                                 0x05,    0x80,0xC8, 0x80,0xC9, 0,10,
                            0x80,0x85,    0x80,0xC8, 0x80,0xC9, 0,10,
        0xF9,0x80,0x00,0x00,0x01,0x80,    1, 0
    }));
}

/*****************************************************************************
 Test the writing of a tad file
 *****************************************************************************/

bool tad_write(std::vector<AlleleDepths> lines, std::string text) {
    std::stringstream buffer;
    Ad adfile("tad:", std::ios_base::out);
    adfile.AddContigs("seq1",100);
    adfile.AddContigs("seq2",1000, "");
    adfile.AddContigs("seq3",10000, "M5:aaaaaaaa");
    adfile.Attach(buffer.rdbuf());
    for(auto line : lines) {
        if(adfile.Write(line) != 0) {
            std::cerr << "  Error: adfile.Write failed\n";
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

bool tad_read(std::string text) {
    std::stringstream buffer;
    Ad adfile("tad:", std::ios_base::in);
    adfile.Attach(buffer.rdbuf());
    adfile.ReadHeader();
}

BOOST_AUTO_TEST_CASE(test_tad_read) {
    BOOST_CHECK(tad_read(
        "@ID\tFF:TAD\tVN:0.1\n"
        "@SQ\tSN:scaffold_1\tLN:100\n"
        "@SQ\tSN:scaffold_2\tLN:200\n"
        "@SQ\tSN:scaffold_3\tLN:300\n"
        "@AD\tID:A\n"
        "@AD\tID:B\n"
        "scaffold_1\t1\tA\t10\t9\n"
        "scaffold_1\t2\tC\t8\t7\n"
        "scaffold_1\t3\tGC\t6,0\t0,2\n"
        "scaffold_2\t1\tT\t4\t0\n"
        "scaffold_3\t1\tA\t0\t4\n"
    ));
}