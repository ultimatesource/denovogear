/*
 * Copyright (c) 2018 Reed A. Cartwright
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
#define BOOST_TEST_MODULE dng::pedigree

#include <dng/pedigree.h>

#include "../testing.h"

#include <dng/detail/rangeio.h>

#include <vector>
#include <boost/optional/optional_io.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>

namespace dng {
template<typename CharType, typename CharTrait>
inline
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const Pedigree::Sex& m) {
    o << static_cast<int>(m);
    return o;
}
} //namespace dng


namespace std {
template<typename CharType, typename CharTrait>
inline
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const std::vector<std::string>& m) {
    o << dng::rangeio::wrap(m);
    return o;
}
} //namespace std


using namespace dng;
using namespace std;
using dng::detail::make_test_range;


BOOST_AUTO_TEST_CASE(test_parse_sex) {
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("."), Pedigree::Sex::Unknown);
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("0"), Pedigree::Sex::Autosomal);
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("1"), Pedigree::Sex::Male);
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("2"), Pedigree::Sex::Female);
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("autosomal"), Pedigree::Sex::Autosomal);
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("male"), Pedigree::Sex::Male);
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("female"), Pedigree::Sex::Female);
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("a"), Pedigree::Sex::Autosomal);
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("m"), Pedigree::Sex::Male);
    BOOST_CHECK_EQUAL(Pedigree::parse_sex("f"), Pedigree::Sex::Female);

    BOOST_CHECK_EQUAL(Pedigree::parse_sex("malewrong"), Pedigree::Sex::Unknown);
}

BOOST_AUTO_TEST_CASE(test_parse_text) {
    using boost::adaptors::filtered;
    using boost::adaptors::transformed;

    const char ped[] = 
        "##PEDNG v1.0\n"
        "A . . 1 =\n"
        "B@founder\t.:0.1\t.\t2\tB1\n"
        "C    A    B    1    C1    C2\n"
        "D A:0.01 B:0.5 2\t=\n"
    ;

    std::vector<std::string> expected_names = {
        "A", "B", "C", "D"
    };
    std::vector<boost::optional<std::string>> expected_dads = {
        boost::none, boost::none, std::string{"A"}, std::string{"A"}
    };
    std::vector<boost::optional<std::string>> expected_moms = {
        boost::none, boost::none, std::string{"B"}, std::string{"B"}
    };    
    std::vector<boost::optional<double>> expected_dad_lens = {
        boost::none, boost::none, boost::none, 0.01
    };
    std::vector<boost::optional<double>> expected_mom_lens = {
        boost::none, boost::none, boost::none, 0.5
    };
    std::vector<std::vector<std::string>> expected_samples = {
        {"A"},{"B1"},{"C1","C2"},{"D"}
    };
    std::vector<std::vector<std::string>> expected_tags = {
        {},{"founder"},{},{}
    };
    auto Male = Pedigree::Sex::Male;
    auto Female = Pedigree::Sex::Female;
    std::vector<Pedigree::Sex> expected_sexes = {
        Male,Female,Male,Female
    };

    Pedigree pedigree;
    BOOST_REQUIRE_NO_THROW(pedigree = Pedigree::parse_text(ped));

    // Check names
    auto test_names = make_test_range(pedigree.table() | 
        transformed(boost::mem_fn(&Pedigree::Member::name)));
    CHECK_EQUAL_RANGES(test_names, expected_names);

    // Check dads
    auto test_dads = make_test_range(pedigree.table() | 
        transformed(boost::mem_fn(&Pedigree::Member::dad)));
    CHECK_EQUAL_RANGES(test_dads, expected_dads);

    // Check moms
    auto test_moms = make_test_range(pedigree.table() | 
        transformed(boost::mem_fn(&Pedigree::Member::mom)));
    CHECK_EQUAL_RANGES(test_moms, expected_moms);

    // Check dad lengths
    auto test_dad_lens = make_test_range(pedigree.table() | 
        transformed(boost::mem_fn(&Pedigree::Member::dad_length)));
    CHECK_EQUAL_RANGES(test_dad_lens, expected_dad_lens);

    // Check mom lengths
    auto test_mom_lens = make_test_range(pedigree.table() | 
        transformed(boost::mem_fn(&Pedigree::Member::mom_length)));
    CHECK_EQUAL_RANGES(test_mom_lens, expected_mom_lens);

    // Check samples
    auto test_samples = make_test_range(pedigree.table() | 
        transformed(boost::mem_fn(&Pedigree::Member::samples)));
    CHECK_EQUAL_RANGES(test_samples, expected_samples);

    // Check tags
    auto test_tags = make_test_range(pedigree.table() | 
        transformed(boost::mem_fn(&Pedigree::Member::tags)));
    CHECK_EQUAL_RANGES(test_tags, expected_tags);

    // Check sex
    auto test_sexes = make_test_range(pedigree.table() | 
        transformed(boost::mem_fn(&Pedigree::Member::sex)));
    CHECK_EQUAL_RANGES(test_sexes, expected_sexes);

    // Check Pedigree Functions
    BOOST_CHECK_EQUAL(pedigree.LookupMemberPosition("B"), 1);
    BOOST_CHECK_EQUAL(pedigree.LookupMember("B")->name, "B");
    BOOST_CHECK_EQUAL(pedigree.GetMember(1).name, "B");
}

BOOST_AUTO_TEST_CASE(test_parse_text_exceptions) {
    Pedigree pedigree;

    const char empty[] = "";
    BOOST_CHECK_THROW(Pedigree::parse_text(empty), std::invalid_argument);

    const char no_pedng[] = "A . . 1 =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(no_pedng), std::invalid_argument);

    const char duplcate_names[] = "##PEDNG v1.0\nA . . 1 =\nA . . 1 =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(duplcate_names), std::invalid_argument);

    const char short_row[] = "##PEDNG v1.0\nA\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(short_row), std::invalid_argument);

    const char no_name[] = "##PEDNG v1.0\n@tag . . 1 =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(no_name), std::invalid_argument);

    const char no_dad_name[] = "##PEDNG v1.0\nA :0.1 C:0.1 1 =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(no_dad_name), std::invalid_argument);

    const char no_dad_length[] = "##PEDNG v1.0\nA B: C:0.1 1 =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(no_dad_length), std::invalid_argument);

    const char bad_dad_length[] = "##PEDNG v1.0\nA B:0.a C:0.1 1 =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(bad_dad_length), std::invalid_argument);

    const char no_mom_name[] = "##PEDNG v1.0\nA B:0.1 :0.1 1 =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(no_mom_name), std::invalid_argument);

    const char no_mom_length[] = "##PEDNG v1.0\nA B:0.1 C: 1 =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(no_mom_length), std::invalid_argument);

    const char bad_mom_length[] = "##PEDNG v1.0\nA B:0.1 C:0.a 1 =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(bad_mom_length), std::invalid_argument);

    const char bad_sex[] = "##PEDNG v1.0\nA B C . =\n";
    BOOST_CHECK_THROW(Pedigree::parse_text(bad_sex), std::invalid_argument);
}
