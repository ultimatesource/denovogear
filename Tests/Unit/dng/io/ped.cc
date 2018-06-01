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
#define BOOST_TEST_MODULE dng::io::ped

#include <dng/io/ped.h>

#include "../../testing.h"

#include <vector>
#include <boost/optional/optional_io.hpp>

using namespace dng;
using namespace dng::io;
using namespace std;

const char simple_ped[] = 
    "##PEDNG v1.0\n"
    "A@founder\t. . \t1    =\t\tA2\n"
;

BOOST_AUTO_TEST_CASE(test_ped) {
    dng::detail::AutoTempFile temp;
    temp.file.write(simple_ped, sizeof(simple_ped)-1);
    temp.file.flush();

    Ped ped;
    BOOST_REQUIRE_NO_THROW(ped = Ped{temp.path.native()});
    
    Pedigree pedigree;
    BOOST_REQUIRE_NO_THROW(pedigree = ped.Parse());

    BOOST_REQUIRE_EQUAL(pedigree.NumberOfMembers(), 1);

    auto const &member = pedigree.GetMember(0);

    BOOST_CHECK_EQUAL(member.name, "A");
    std::vector<string> expected_tags = {"founder"};
    CHECK_EQUAL_RANGES(member.tags, expected_tags);
    BOOST_CHECK_EQUAL(member.dad, boost::none);
    BOOST_CHECK_EQUAL(member.dad_length, boost::none);
    BOOST_CHECK_EQUAL(member.mom, boost::none);
    BOOST_CHECK_EQUAL(member.mom_length, boost::none);
    BOOST_CHECK_EQUAL(member.sex, Pedigree::Sex::Male);
    std::vector<string> expected_samples = {"A", "A2"};
    CHECK_EQUAL_RANGES(member.samples, expected_samples);
}

BOOST_AUTO_TEST_CASE(test_parse_ped) {
    dng::detail::AutoTempFile temp;
    temp.file.write(simple_ped, sizeof(simple_ped)-1);
    temp.file.flush();
    
    Pedigree pedigree;
    BOOST_REQUIRE_NO_THROW(pedigree = parse_ped(temp.path.native()));

    BOOST_REQUIRE_EQUAL(pedigree.NumberOfMembers(), 1);

    auto const &member = pedigree.GetMember(0);

    BOOST_CHECK_EQUAL(member.name, "A");
    std::vector<string> expected_tags = {"founder"};
    CHECK_EQUAL_RANGES(member.tags, expected_tags);
    BOOST_CHECK_EQUAL(member.dad, boost::none);
    BOOST_CHECK_EQUAL(member.dad_length, boost::none);
    BOOST_CHECK_EQUAL(member.mom, boost::none);
    BOOST_CHECK_EQUAL(member.mom_length, boost::none);
    BOOST_CHECK_EQUAL(member.sex, Pedigree::Sex::Male);
    std::vector<string> expected_samples = {"A", "A2"};
    CHECK_EQUAL_RANGES(member.samples, expected_samples);
}
