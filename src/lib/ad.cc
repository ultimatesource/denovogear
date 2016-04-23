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
#include <dng/utility.h>
#include <dng/io/ad.h>
#include <dng/io/utility.h>

#include <algorithm>
#include <cstdlib>

#include <endian.h>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/std_pair.hpp>

using namespace dng::pileup;

const AlleleDepths::type_info_t
    AlleleDepths::type_info_table[128] = {
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

constexpr int type_info_table_length = 128;
static_assert(sizeof(AlleleDepths::type_info_table) / sizeof(AlleleDepths::type_info_t) == type_info_table_length,
    "AlleleDepths::type_info_table does not have 128 elements.");

struct tad_label_lookup_t {
    boost::spirit::qi::tst<char, int> tst_;
    tad_label_lookup_t() {
        for(int i=0; i<type_info_table_length; ++i) {
            auto & slot = AlleleDepths::type_info_table[i];
            tst_.add(slot.label_upper, slot.label_upper+strlen(slot.label_upper),i);
        }
    }
    int operator()(const std::string &s) {
        auto first = s.cbegin();
        auto last = s.cend();

        int *p = tst_.find(first,last, [](char ch) -> char { return std::toupper(ch); });
        if(first != last || p == nullptr) {
            return -1;
        }
        return *p;
    }
};

int8_t AlleleDepths::LookupType(const std::string& ss) {
    static tad_label_lookup_t db;
    return db(ss);
}

int8_t AlleleDepths::LookupType(const std::vector<std::size_t> &indexes, bool ref_is_n) {
    static size_t starts[] = {0, 4, 16, 40, 64};

    assert(1 <= indexes.size() && indexes.size() <= 4);
    auto start = &type_info_table[0] + starts[indexes.size()-1] + (ref_is_n ? 64 : 0);
    auto end =   &type_info_table[0] + starts[indexes.size()] + (ref_is_n ? 64 : 0);

    auto it = std::find_if(start, end, [&](const type_info_t& type) {
        for(size_t u = 0; u < indexes.size(); ++u) {
            if(type.indexes[u] != indexes[u]) {
                return false;
            }
        }
        return true;
    });
    return (it != end) ? it->type : -1;
}

int ntf8_put32(uint32_t n, char *out, size_t count);
int ntf8_put64(uint64_t n, char *out, size_t count);
int ntf8_get32(const char *in, size_t count, uint32_t *r);
int ntf8_get64(const char *in, size_t count, uint64_t *r);

void dng::io::Ad::Clear() {
    contigs_.clear();
    id_.name.clear();
    id_.version = 0;
    id_.attributes.clear();
    extra_headers_.clear();
    contig_tst_.clear();    
}

void tab_append(std::string *out, const std::string &in, char d='\t') {
    assert(out != nullptr);
    if(out->empty()) {
        *out = in;
    } else {
        *out += d;
        *out += in;
    }
}

int dng::io::Ad::ReadHeaderTad() {
    using namespace boost;
    using namespace std;
    using boost::algorithm::starts_with;
    using boost::algorithm::iequals;
    using boost::algorithm::to_upper;
    namespace ss = boost::spirit::standard;
    namespace qi = boost::spirit::qi;
    using ss::space;
    using qi::uint_;
    using qi::lit;
    using qi::phrase_parse;
    qi::uint_parser<uint8_t> uint8_;

    // Construct the tokenizer
    typedef tokenizer<char_separator<char>,
            std::istreambuf_iterator<char>> tokenizer;
    char_separator<char> sep("\t", "\n");
    tokenizer tokens(istreambuf_range(stream_), sep);

    // Store all the tokens in the header until we reach the first data line
    vector<string> store;
    for(auto it = tokens.begin(); it != tokens.end(); ++it) {
        store.push_back(*it);
        if(*it == "\n" && stream_.peek() != '@') {
            break;
        }
    }
    // Clear header information
    Clear();

    // Setup header information
    auto it = store.begin();
    if(it == store.end() || *it != "@ID") {
        throw(runtime_error("Unable to parse TAD header: @ID missing from first line."));
    }
    // Parse @ID line which identifies the file format and version.
    for(++it; it != store.end(); ++it) {
        if(*it == "\n") {
            ++it;
            break;
        } else if(starts_with(*it, "FF:")) {
            id_.name = it->substr(3);
            to_upper(id_.name);
        } else if(starts_with(*it, "VN:")) {
            const char *first = it->c_str()+3;
            const char *last  = it->c_str()+it->length();
            std:pair<uint8_t,uint8_t> bytes;
            if(!phrase_parse(first, last, uint8_ >> ( '.' >> uint8_ || lit('\0')), space, bytes) || first != last) {
                throw(runtime_error("Unable to parse TAD header: File format version information '" + *it + "' not major.minor."));
            }
            id_.version = ((uint16_t)bytes.first << 8) | bytes.second;
        } else {
            tab_append(&id_.attributes, *it);
        }
    }
    if(id_.name != "TAD") {
        throw(runtime_error("Unable to parse TAD header: Unknown file format '" + id_.name + "'"));
    }
    for(; it != store.end(); ++it) {
        if(*it == "@SQ") {
            contig_t sq;
            for(++it; it != store.end() && *it != "\n"; ++it) {
                if(starts_with(*it, "SN:")) {
                    sq.name = it->substr(3);
                } else if(starts_with(*it, "LN:")) {
                    const char *first = it->c_str()+3;
                    const char *last  = it->c_str()+it->length();

                    if(!phrase_parse(first, last, uint_ , space, sq.length) || first != last) {
                        throw(runtime_error("Unable to parse TAD header: Unable to parse sequence length from header attribute '" + *it + "'"));
                    }
                } else {
                    tab_append(&sq.attributes, *it);
                }
            }
            contigs_.push_back(std::move(sq));
        } else if(*it == "@AD") {
            library_t ad;
            for(++it; it != store.end() && *it != "\n"; ++it) {
                if(starts_with(*it, "ID:")) {
                    ad.name = it->substr(3);
                } else if(starts_with(*it, "SM:")) {
                    ad.sample = it->substr(3);
                } else {
                    tab_append(&ad.attributes, *it);
                }
            }
            libraries_.push_back(std::move(ad));
        } else if(*it == "@ID") {
            throw(runtime_error("Unable to parse TAD header: @ID can only be specified once."));
        } else if((*it)[0] == '@') {
            // store all the information for this line into the extra_headers_ field
            extra_headers_.push_back(*it);
            for(++it; it != store.end() && *it != "\n"; ++it) {
                extra_headers_.push_back(*it);
            }
            extra_headers_.push_back("\n");
        } else {
            throw(runtime_error("Unable to parse TAD header: Something went wrong. "
                "Expected the @TAG at the beginning of a header line but got '" + *it +"'."));            
        }
    }
    return 0;
}

int dng::io::Ad::ReadTad(AlleleDepths *pline) {
    using namespace std;
    using namespace boost;
    namespace ss = boost::spirit::standard;
    namespace qi = boost::spirit::qi;
    using ss::space;
    using qi::uint_;
    using qi::phrase_parse;

    if(pline == nullptr) {
        // discard this line
        while(stream_ && stream_.get() != '\n') {
            /*noop*/;
        }
        return 0;
    }
    // Construct the tokenizer
    typedef tokenizer<char_separator<char>,
            std::istreambuf_iterator<char>> tokenizer;
    char_separator<char> sep("\t", "\n");
    tokenizer tokens(istreambuf_range(stream_), sep);
    vector<string> store;
    for(auto it = tokens.begin(); it != tokens.end(); ++it) {
        if(*it == "\n") {
            break;
        }        
        store.push_back(*it);
    }
    // check to make sure store has the proper size
    if(store.size() != 3+num_libraries_) {
        throw(runtime_error("Unable to parse TAD record: Expected " + utility::to_pretty(3+num_libraries_) 
            + " columns. Found " + utility::to_pretty(store.size()) + " columns."));
    }

    // parse contig
    auto first = store[0].cbegin();
    auto last = store[0].cend();

    int *p = contig_tst_.find(first,last);
    if(first != last || p == nullptr) {
        throw(runtime_error("Unable to parse TAD record: Contig Name '" + store[0] + "' is unknown."));
    }
    int contig_num = *p;
    int position;
    // parse position
    first = store[1].cbegin();
    last  = store[1].cend();

    if(!phrase_parse(first, last, uint_, space, position) || first != last) {
        throw(runtime_error("Unable to parse TAD record: Unable to parse position '" + store[1] + "' as integer."));
    }
    pline->location(utility::make_location(contig_num,position-1));

    // parse type
    int8_t ty = AlleleDepths::LookupType(store[2]);
    if(ty == -1) {
        throw(runtime_error("Unable to parse TAD record: Type '" + store[2] + "' is unknown."));        
    }
    pline->resize(ty,num_libraries_);

    // parse data
    vector<int> data;
    bool success;
    for(size_t i=3;i<store.size();++i) {
        tie(data,success) = utility::parse_int_list(store[i]);
        for(size_t j=0;j<data.size();++j) {
            (*pline)(j,i-3) = data[j];
        }
    }
    return 0;
}

int dng::io::Ad::WriteTad(const AlleleDepths& line) {
    using dng::utility::location_to_contig;
    using dng::utility::location_to_position;
    typedef AlleleDepths::size_type size_type;

    if(line.data_size() == 0) {
        return 0;
    }
    const location_t loc = line.location();
    const size_type nlib = line.num_libraries();
    const size_type nnuc = line.num_nucleotides();
    auto contig_num = location_to_contig(loc);
    if(contig_num >= contigs_.size()) {
        return 1;
    }
    stream_ << contigs_[contig_num].name << '\t';
    stream_ << location_to_position(loc)+1 << '\t';
    stream_ << line.type_info().label_upper;
    for(AlleleDepths::size_type y = 0; y != nlib; ++y) {
        stream_  << '\t';
        size_type x = 0;
        for(; x != nnuc-1; ++x) {
            stream_  << line(x,y) << ',';
        }
        stream_  << line(x,y);
    }
    stream_ << '\n';
    return 0;
}

int dng::io::Ad::ReadHeaderAd() {
    assert(false); // not implemented yet
    return 1;
}

int dng::io::Ad::ReadAd(AlleleDepths *pline) {
    assert(false); // not implemented yet
    return 1;
}

int dng::io::Ad::WriteAd(const AlleleDepths& line) {
    static_assert(sizeof(AlleleDepths::type_info_table) / sizeof(AlleleDepths::type_info_t) == 128,
        "AlleleDepths::type_info_table does not have 128 elements.");
    static_assert(sizeof(location_t) == 8, "location_t is not a 64-bit value.");

    using dng::utility::location_to_contig;
    using dng::utility::location_to_position;
    typedef AlleleDepths::size_type size_type;

    const size_type nlib = line.num_libraries();
    if(line.data_size() == 0) {
        return 0;
    }
    location_t loc = line.location();
    // The Ad format only holds up to 2^25-2 contigs 
    if(location_to_contig(loc) > 0x01FFFFFE) {
        return 1;
    }

    if(counter_ == 0 || location_to_contig(loc) != location_to_contig(last_location_)) {
        // use absolute positioning
        // Adjust the contig number by 1
        last_location_ = loc;
        loc += (1LL << 32);
    } else {
        // Make sure we are sorted here
        assert(loc > last_location_);
        location_t diff = loc - last_location_ - 1;
        last_location_ = loc;
        loc = diff;
    }
    // set counter
    counter_ += 1;

    // Fetch type and check if it is positive
    int8_t type = line.type();
    assert(type >= 0);
    int64_t u = (loc << 7) | type; 
    char buffer[9];
    // write out contig, position, and location
    int sz = ntf8_put64(u, buffer,sizeof(buffer));
    stream_.write(buffer, sz);
    // write out data
    for(size_type i = 0; i < line.data_size(); ++i) {
        sz = ntf8_put32(line.data()[i], buffer,sizeof(buffer));
        stream_.write(buffer,sz);
    }

    return 0;
}

/* NTF8 Format
NTF8 is a multi-byte, big-endian numeric format.
The binary format consists of a prefix of 0 or more 1s,
followed by a zero, followed by the value in big-endian format.
The maximum amount of data that a single NTF8 number can hold
is 64-bits. The number of leading 1s in the first byte tells
the user how many extra bytes need to be read to decode
the value.
*/

int ntf8_put32(uint32_t n, char *out, size_t count) {
    if(n <= 0x7F) {
        // 0bbb bbbb
        if(count < 1) {
            return -1;
        }
        *out = n;
        return 1;
    } else if(n <= 0x3FFF) {
        // 10bb bbbb bbbb bbbb
        if(count < 2) {
            return -1;
        }
    	uint16_t u = htobe16(n | 0x8000);
    	memcpy(out, &u, 2);
        return 2;
    } else if(n <= 0x1FFFFF) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        if(count < 3) {
            return -1;
        }
        uint32_t u = htobe32((n | 0xC00000) << 8);
    	memcpy(out, &u, 3);
        return 3;
    } else if(n <= 0x0FFFFFFF) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 4) {
            return -1;
        }
        uint32_t u = htobe32(n | 0xE0000000);
    	memcpy(out, &u, 4);
        return 4;
    } else {
        // 1111 0000 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 5) {
            return -1;
        }
        *out = 0xF0;
        uint32_t u = htobe32(n);
        memcpy(out+1, &u, 4);
        return 5;
    }
}

int ntf8_get32(const char *in, size_t count, uint32_t *r) {
	assert(r != nullptr);
    if(count == 0) {
        return -1;
    }
    uint8_t x = in[0];
    if(x < 0x80) {
        // 0bbb bbbb     
        *r = x;
        return 1;
    } else if(x < 0xC0) {
        // 10bb bbbb bbbb bbbb
        if(count < 2) {
            return -1;
        }
        uint16_t u = 0;
        memcpy(&u,in,2);
        *r = be16toh(u) & 0x3FFF;
        return 2;
    } else if(x < 0xE0) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        if(count < 3) {
            return -1;
        }
        uint32_t u = 0;
        memcpy(&u,in,3);
        *r = (be32toh(u) >> 8) & 0x1FFFFF;
        return 3;
    } else if(x < 0xF0) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 4) {
            return -1;
        }
        uint32_t u = 0;
        memcpy(&u,in,4);
        *r = be32toh(u) & 0x0FFFFFFF;
        return 4;
    } else {
        assert(x == 0xF0); // If this fails then we may be reading a 64-bit number
        // 1111 0000 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 5) {
            return -1;
        }        
        uint32_t u = 0;
        memcpy(&u,in+1,4);
        *r = be32toh(u);
        return 5;
    }
}

int ntf8_put64(uint64_t n, char *out, size_t count) {
    if(n <= 0x7F) {
        // 0bbb bbbb
        if(count < 1) {
            return -1;
        }
        *out = n;
        return 1;
    } else if(n <= 0x3FFF) {
        // 10bb bbbb bbbb bbbb
        if(count < 2) {
            return -1;
        }
    	uint16_t u = htobe16(n | 0x8000);
    	memcpy(out, &u, 2);
        return 2;
    } else if(n <= 0x1FFFFF) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        if(count < 3) {
            return -1;
        }
        uint32_t u = htobe32((n | 0xC00000) << 8);
    	memcpy(out, &u, 3);
        return 3;
    } else if(n <= 0x0FFFFFFF) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 4) {
            return -1;
        }
        uint32_t u = htobe32(n | 0xE0000000);
    	memcpy(out, &u, 4);
        return 4;
    } else if(n <= 0x07FFFFFFFF) {
        // 1111 0bbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 5) {
            return -1;
        }
        uint64_t u = htobe64((n | 0xF000000000) << 24);
        memcpy(out, &u, 5);
        return 5;
    } else if(n <= 0x03FFFFFFFFFF) {
        // 1111 10bb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 6) {
            return -1;
        }
        uint64_t u = htobe64((n | 0xF80000000000) << 16);
        memcpy(out, &u, 6);
        return 6;
    } else if(n <= 0x01FFFFFFFFFFFF) {
        // 1111 110b bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 7) {
            return -1;
        }
        uint64_t u = htobe64((n | 0xFC000000000000) << 8);
        memcpy(out, &u, 7);
        return 7;
    } else if(n <= 0x00FFFFFFFFFFFFFF) {
        // 1111 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 8) {
            return -1;
        }
        uint64_t u = htobe64( n | 0xFE00000000000000);
        memcpy(out, &u, 8);
        return 8;
    } else {
        // 1111 1111 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 9) {
            return -1;
        }
        *out = 0xFF;
        uint64_t u = htobe64(n);
        memcpy(out+1, &u, 8);
        return 9;
    }
}

int ntf8_get64(const char *in, size_t count, uint64_t *r) {
    assert(r != nullptr);
    if(count == 0) {
        return -1;
    }
    uint8_t x = in[0];
    if(x < 0x80) {
        // 0bbb bbbb
        *r = x;
        return 1;
    } else if(x < 0xC0) {
        // 10bb bbbb bbbb bbbb
        if(count < 2) {
            return -1;
        }
        uint16_t u = 0;
        memcpy(&u,in,2);
        *r = be16toh(u) & 0x3FFF;
        return 2;
    } else if(x < 0xE0) {
        // 110b bbbb bbbb bbbb bbbb bbbb
        if(count < 3) {
            return -1;
        }
        uint32_t u = 0;
        memcpy(&u,in,3);
        *r = (be32toh(u) >> 8) & 0x1FFFFF;
        return 3;
    } else if(x < 0xF0) {
        // 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 4) {
            return -1;
        }
        uint32_t u = 0;
        memcpy(&u,in,4);
        *r = be32toh(u) & 0x0FFFFFFF;
        return 4;
    } else if(x < 0xF8) {
        // 1111 0bbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 5) {
            return -1;
        }
        uint64_t u = 0;
        memcpy(&u,in,5);
        *r = (be64toh(u) >> 24) & 0x07FFFFFFFF;
        return 5;
    } else if(x < 0xFC) {
        // 1111 10bb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 6) {
            return -1;
        }
        uint64_t u = 0;
        memcpy(&u,in,6);
        *r = (be64toh(u) >> 16) & 0x03FFFFFFFFFF;
        return 6;
    } else if(x < 0xFE) {
        // 1111 110b bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 7) {
            return -1;
        }
        uint64_t u = 0;
        memcpy(&u,in,7);
        *r = (be64toh(u) >> 8) & 0x01FFFFFFFFFFFF;
        return 7;
    } else if(x < 0xFF) {
        // 1111 1110 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 8) {
            return -1;
        }
        uint64_t u = 0;
        memcpy(&u,in,8);
        *r = be64toh(u) & 0x00FFFFFFFFFFFFFF;
        return 8;
    } else {
        // 1111 1111 bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb bbbb
        if(count < 9) {
            return -1;
        }
        uint64_t u = 0;
        memcpy(&u,in+1,8);
        *r = be64toh(u);
        return 9;
    }
}
