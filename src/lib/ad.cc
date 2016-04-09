/*
 * Copyright (c) 2015-2016 Reed A. Cartwright
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

#include <dng/depths.h>

#include <algorithm>
#include <endian.h>

const AlleleDepths::type_info_t AlleleDepths::type_info_table[128] = {
    {0,   1, "A",     "a",     0, {0,3,2,1}},
    {1,   1, "C",     "c",     1, {1,3,2,0}},
    {2,   1, "G",     "g",     2, {2,3,1,0}},
    {3,   1, "T",     "t",     3, {3,2,1,0}},
    {4,   2, "AC",    "ac",    0, {0,1,3,2}},
    {5,   2, "AG",    "ag",    0, {0,2,3,1}},
    {6,   2, "AT",    "at",    0, {0,3,2,1}},
    {7,   2, "CA",    "ca",    1, {1,0,3,2}},
    {8,   2, "CG",    "cg",    1, {1,2,3,0}},
    {9,   2, "CT",    "ct",    1, {1,3,2,0}},
    {10,  2, "GA",    "ga",    2, {2,0,3,1}},
    {11,  2, "GC",    "gc",    2, {2,1,3,0}},
    {12,  2, "GT",    "gt",    2, {2,3,1,0}},
    {13,  2, "TA",    "ta",    3, {3,0,2,1}},
    {14,  2, "TC",    "tc",    3, {3,1,2,0}},
    {15,  2, "TG",    "tg",    3, {3,2,1,0}},
    {16,  3, "ACG",   "acg",   0, {0,1,2,3}},
    {17,  3, "ACT",   "act",   0, {0,1,3,2}},
    {18,  3, "AGC",   "agc",   0, {0,2,1,3}},
    {19,  3, "AGT",   "agt",   0, {0,2,3,1}},
    {20,  3, "ATC",   "atc",   0, {0,3,1,2}},
    {21,  3, "ATG",   "atg",   0, {0,3,2,1}},
    {22,  3, "CAG",   "cag",   1, {1,0,2,3}},
    {23,  3, "CAT",   "cat",   1, {1,0,3,2}},
    {24,  3, "CGA",   "cga",   1, {1,2,0,3}},
    {25,  3, "CGT",   "cgt",   1, {1,2,3,0}},
    {26,  3, "CTA",   "cta",   1, {1,3,0,2}},
    {27,  3, "CTG",   "ctg",   1, {1,3,2,0}},
    {28,  3, "GAC",   "gac",   2, {2,0,1,3}},
    {29,  3, "GAT",   "gat",   2, {2,0,3,1}},
    {30,  3, "GCA",   "gca",   2, {2,1,0,3}},
    {31,  3, "GCT",   "gct",   2, {2,1,3,0}},
    {32,  3, "GTA",   "gta",   2, {2,3,0,1}},
    {33,  3, "GTC",   "gtc",   2, {2,3,1,0}},
    {34,  3, "TAC",   "tac",   3, {3,0,1,2}},
    {35,  3, "TAG",   "tag",   3, {3,0,2,1}},
    {36,  3, "TCA",   "tca",   3, {3,1,0,2}},
    {37,  3, "TCG",   "tcg",   3, {3,1,2,0}},
    {38,  3, "TGA",   "tga",   3, {3,2,0,1}},
    {39,  3, "TGC",   "tgc",   3, {3,2,1,0}},
    {40,  4, "ACGT",  "acgt",  0, {0,1,2,3}},
    {41,  4, "ACTG",  "actg",  0, {0,1,3,2}},
    {42,  4, "AGCT",  "agct",  0, {0,2,1,3}},
    {43,  4, "AGTC",  "agtc",  0, {0,2,3,1}},
    {44,  4, "ATCG",  "atcg",  0, {0,3,1,2}},
    {45,  4, "ATGC",  "atgc",  0, {0,3,2,1}},
    {46,  4, "CAGT",  "cagt",  1, {1,0,2,3}},
    {47,  4, "CATG",  "catg",  1, {1,0,3,2}},
    {48,  4, "CGAT",  "cgat",  1, {1,2,0,3}},
    {49,  4, "CGTA",  "cgta",  1, {1,2,3,0}},
    {50,  4, "CTAG",  "ctag",  1, {1,3,0,2}},
    {51,  4, "CTGA",  "ctga",  1, {1,3,2,0}},
    {52,  4, "GACT",  "gact",  2, {2,0,1,3}},
    {53,  4, "GATC",  "gatc",  2, {2,0,3,1}},
    {54,  4, "GCAT",  "gcat",  2, {2,1,0,3}},
    {55,  4, "GCTA",  "gcta",  2, {2,1,3,0}},
    {56,  4, "GTAC",  "gtac",  2, {2,3,0,1}},
    {57,  4, "GTCA",  "gtca",  2, {2,3,1,0}},
    {58,  4, "TACG",  "tacg",  3, {3,0,1,2}},
    {59,  4, "TAGC",  "tagc",  3, {3,0,2,1}},
    {60,  4, "TCAG",  "tcag",  3, {3,1,0,2}},
    {61,  4, "TCGA",  "tcga",  3, {3,1,2,0}},
    {62,  4, "TGAC",  "tgac",  3, {3,2,0,1}},
    {63,  4, "TGCA",  "tgca",  3, {3,2,1,0}},
    {64,  1, "NA",    "na",    4, {0,3,2,1}},
    {65,  1, "NC",    "nc",    4, {1,3,2,0}},
    {66,  1, "NG",    "ng",    4, {2,3,1,0}},
    {67,  1, "NT",    "nt",    4, {3,2,1,0}},
    {68,  2, "NAC",   "nac",   4, {0,1,3,2}},
    {69,  2, "NAG",   "nag",   4, {0,2,3,1}},
    {70,  2, "NAT",   "nat",   4, {0,3,2,1}},
    {71,  2, "NCA",   "nca",   4, {1,0,3,2}},
    {72,  2, "NCG",   "ncg",   4, {1,2,3,0}},
    {73,  2, "NCT",   "nct",   4, {1,3,2,0}},
    {74,  2, "NGA",   "nga",   4, {2,0,3,1}},
    {75,  2, "NGC",   "ngc",   4, {2,1,3,0}},
    {76,  2, "NGT",   "ngt",   4, {2,3,1,0}},
    {77,  2, "NTA",   "nta",   4, {3,0,2,1}},
    {78,  2, "NTC",   "ntc",   4, {3,1,2,0}},
    {79,  2, "NTG",   "ntg",   4, {3,2,1,0}},
    {80,  3, "NACG",  "nacg",  4, {0,1,2,3}},
    {81,  3, "NACT",  "nact",  4, {0,1,3,2}},
    {82,  3, "NAGC",  "nagc",  4, {0,2,1,3}},
    {83,  3, "NAGT",  "nagt",  4, {0,2,3,1}},
    {84,  3, "NATC",  "natc",  4, {0,3,1,2}},
    {85,  3, "NATG",  "natg",  4, {0,3,2,1}},
    {86,  3, "NCAG",  "ncag",  4, {1,0,2,3}},
    {87,  3, "NCAT",  "ncat",  4, {1,0,3,2}},
    {88,  3, "NCGA",  "ncga",  4, {1,2,0,3}},
    {89,  3, "NCGT",  "ncgt",  4, {1,2,3,0}},
    {90,  3, "NCTA",  "ncta",  4, {1,3,0,2}},
    {91,  3, "NCTG",  "nctg",  4, {1,3,2,0}},
    {92,  3, "NGAC",  "ngac",  4, {2,0,1,3}},
    {93,  3, "NGAT",  "ngat",  4, {2,0,3,1}},
    {94,  3, "NGCA",  "ngca",  4, {2,1,0,3}},
    {95,  3, "NGCT",  "ngct",  4, {2,1,3,0}},
    {96,  3, "NGTA",  "ngta",  4, {2,3,0,1}},
    {97,  3, "NGTC",  "ngtc",  4, {2,3,1,0}},
    {98,  3, "NTAC",  "ntac",  4, {3,0,1,2}},
    {99,  3, "NTAG",  "ntag",  4, {3,0,2,1}},
    {100, 3, "NTCA",  "ntca",  4, {3,1,0,2}},
    {101, 3, "NTCG",  "ntcg",  4, {3,1,2,0}},
    {102, 3, "NTGA",  "ntga",  4, {3,2,0,1}},
    {103, 3, "NTGC",  "ntgc",  4, {3,2,1,0}},
    {104, 4, "NACGT", "nacgt", 4, {0,1,2,3}},
    {105, 4, "NACTG", "nactg", 4, {0,1,3,2}},
    {106, 4, "NAGCT", "nagct", 4, {0,2,1,3}},
    {107, 4, "NAGTC", "nagtc", 4, {0,2,3,1}},
    {108, 4, "NATCG", "natcg", 4, {0,3,1,2}},
    {109, 4, "NATGC", "natgc", 4, {0,3,2,1}},
    {110, 4, "NCAGT", "ncagt", 4, {1,0,2,3}},
    {111, 4, "NCATG", "ncatg", 4, {1,0,3,2}},
    {112, 4, "NCGAT", "ncgat", 4, {1,2,0,3}},
    {113, 4, "NCGTA", "ncgta", 4, {1,2,3,0}},
    {114, 4, "NCTAG", "nctag", 4, {1,3,0,2}},
    {115, 4, "NCTGA", "nctga", 4, {1,3,2,0}},
    {116, 4, "NGACT", "ngact", 4, {2,0,1,3}},
    {117, 4, "NGATC", "ngatc", 4, {2,0,3,1}},
    {118, 4, "NGCAT", "ngcat", 4, {2,1,0,3}},
    {119, 4, "NGCTA", "ngcta", 4, {2,1,3,0}},
    {120, 4, "NGTAC", "ngtac", 4, {2,3,0,1}},
    {121, 4, "NGTCA", "ngtca", 4, {2,3,1,0}},
    {122, 4, "NTACG", "ntacg", 4, {3,0,1,2}},
    {123, 4, "NTAGC", "ntagc", 4, {3,0,2,1}},
    {124, 4, "NTCAG", "ntcag", 4, {3,1,0,2}},
    {125, 4, "NTCGA", "ntcga", 4, {3,1,2,0}},
    {126, 4, "NTGAC", "ntgac", 4, {3,2,0,1}},
    {127, 4, "NTGCA", "ntgca", 4, {3,2,1,0}}
};

static int8_t AlleleDepth::LookupType(const std::vector<std::size_t> &indexes, bool ref_is_n) {
    static size_t starts[] = {0, 4, 16, 40, 64};

    assert(1 <= indexes.size() && indexes.size() <= 4);
    size_t start = type_info_table + starts[indexes.size()-1] + (ref_is_n ? 64 : 0);
    size_t end =   type_info_table + starts[indexes.size()] + (ref_is_n ? 64 : 0);

    auto it = std::find_if(start, end, [&](const type_info_t& type) {
        for(size_t u = 0; u < indexes.size(); ++u) {
            if(type.indexes[u] != indexes[u]) {
                return false;
            }
        }
        return true;
    });
    return (it != end) it->type : -1;
}


int itf8_put(char *cp, int32_t val);
int ltf8_put(char *cp, int64_t val);
int ad_write_line(dng::location_t last_location, const dng::io::Ad::AlleleDepths& line);
int tad_write_line(const dng::io::Ad::AlleleDepths& line);

int dng::io::Ad::Write(const AlleleDepths& line) {
    if(is_binary_ad_) {
        return ad_write_line(last_location_, line);
    } else {
        return tad_write_line(line);
    }
}

int ntf8_put32(char *out, uint32_t n) {
    if(n <= 0x7F) {
        // 0bbb bbbb
        *out = n;
        return 1;
    } else if(n <= 0x3FFF) {
        // 10bb bbbb bbbb bbbb
    	uint16_t u = htobe16(n | 0x8000);
    	memcpy(out, &u, 2);
        return 2;
    } else if(n <= 0x1FFFFF) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        uint32_t u = htobe32((n << 8) | 0xC0000000);
    	memcpy(out, &u, 3);
        return 3;
    } else if(n <= 0x0FFFFFFF) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint32_t u = htobe32(n | 0xE0000000);
    	memcpy(out, &u, 4);
        return 4;
    } else {
        // 1111 0000 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        *out = 0xF0;
        uint32_t u = htobe32(n);
        memcpy(out+1, &u, 4);
        return 5;
    }
}

int ntf8_get32(char *in, uint32_t *r) {
	assert(r != nullptr)
    uint8_t x = in[0];
    if(x < 0x80) {
        // 0bbb bbbb
        *r = x;
        return 1;
    } else if(x < 0xC0) {
        // 10bb bbbb bbbb bbbb
        uint16_t u = 0;
        memcpy(&u,in,2);
        *r = be16toh(u) & 0x3FFF;
        return 2;
    } else if(x < 0xE0) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        uint32_t u = 0;
        memcpy(&u,in,3);
        *r = (be32toh(u) & 0x1FFFFFFF) >> 8;
        return 3;
    } else if(x < 0xF0) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint32_t u = 0;
        memcpy(&u,in,4);
        *r = be32toh(u) & 0x0FFFFFFF;
        return 4;
    } else {
        // 1111 0000 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint32_t u = 0;
        memcpy(&u,in+1,4);
        return 5;
    }
}

int ntf8_put64(char *out, uint64_t n) {
    if(n <= 0x7F) {
        // 0bbb bbbb
        *out = n;
        return 1;
    } else if(n <= 0x3FFF) {
        // 10bb bbbb bbbb bbbb
    	uint16_t u = htobe16(n | 0x8000);
    	memcpy(out, &u, 2);
        return 2;
    } else if(n <= 0x1FFFFF) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        uint32_t u = htobe32((n << 8) | 0xC0000000);
    	memcpy(out, &u, 3);
        return 3;
    } else if(n <= 0x0FFFFFFF) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint32_t u = htobe32(n | 0xE0000000);
    	memcpy(out, &u, 4);
        return 4;
    } else if(n <= 0x07FFFFFFFF) {
        // 1111 0bbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t u = htobe64((n << 24) | 0xF000000000000000);
        memcpy(out, &u, 5);
        return 5;
    } else if(n <= 0x03FFFFFFFFFF) {
        // 1111 10bb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t u = htobe64((n << 16) | 0xF800000000000000);
        memcpy(out, &u, 6);
        return 6;
    } else if(n <= 0x01FFFFFFFFFFFF) {
        // 1111 110b bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t u = htobe64((n << 8) | 0xFC00000000000000);
        memcpy(out, &u, 7);
        return 7;
    } else if(n <= 0x00FFFFFFFFFFFFFF) {
        // 1111 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t u = htobe64( n | 0xFE00000000000000);
        memcpy(out, &u, 8);
        return 8;
    } else {
        // 1111 1111 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        *out = 0xFF;
        uint64_t u = htobe64(n);
        memcpy(out+1, &u, 8);
        return 9;
    }
}

std::pair<uint64_t,int> ntf8_get64(char *in) {
    char x = in[0];
    if(x < 0x80) {
        // 0bbb bbbb
        uint64_t r = x;
        return {r,1};
    } else if(x < 0xC0) {
        // 10bb bbbb bbbb bbbb
        uint64_t r = ((x << 8) | in[1]) & 0x03FFFF;
        return {r,2};
    } else if(x < 0xE0) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        uint64_t r = ((x << 16) | (in[1] << 8) | in[2]) & 0x01FFFFFF;
        return {r,3};
    } else if(x < 0xF0) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t r = ((x << 24) | (in[1] << 16) | (in[2] << 8) | in[3]) & 0x0FFFFFFF;
        return {r,4};
    } else if(x < 0xF8) {
        // 1111 0bbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t r = ((x << 32) | (in[1] << 24) | (in[2] << 16) | (in[3] << 8) | in[4]) & 0x07FFFFFFFF;
        return {r,5};
    } else if(x < 0xFC) {
        // 1111 10bb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t r = x << 40;
        r |= in[1] << 32;
        r |= in[2] << 24;
        r |= in[3] << 16;
        r |= in[4] << 8;
        r |= in[5];
        r &= 0x03FFFFFFFFFF;
        return {r,6};
    } else if(x < 0xFE) {
        // 1111 110b bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t r = x << 48;
        r |= in[1] << 40;
        r |= in[2] << 32;
        r |= in[3] << 24;
        r |= in[4] << 16;
        r |= in[5] << 8;
        r |= in[6];
        r &= 0x01FFFFFFFFFFFF;
        return {r,7};
    } else if(x < 0xFF) {
        // 1111 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t r = x << 56;
        r |= in[1] << 48;
        r |= in[2] << 40;
        r |= in[3] << 32;
        r |= in[4] << 24;
        r |= in[5] << 16;
        r |= in[6] << 8;
        r |= in[7];
        r &= 0x00FFFFFFFFFFFFFF;
        return {r,8};
    } else {
        // 1111 1111 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        uint64_t r = in[1] << 56;
        r |= in[2] << 48;
        r |= in[3] << 40;
        r |= in[4] << 32;
        r |= in[5] << 24;
        r |= in[6] << 16;
        r |= in[7] << 8;
        r |= in[8];
       return {r,9};
    }
}
