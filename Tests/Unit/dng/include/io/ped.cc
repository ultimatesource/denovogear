/*
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
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


#define BOOST_TEST_MODULE dng::include::io::ped

#include <dng/io/ped.h>

#include <iostream>
#include <vector>

#include "../../..//boost_test_helper.h"
#include "../../fixture_trio_workspace.h"

namespace dng {
namespace io {
//BOOST_FIXTURE_TEST_CASE(test_constructor, TrioWorkspace) {


void boost_check_equal_member(Pedigree::Member expected, Pedigree::Member actual){

    BOOST_CHECK_EQUAL(expected.fam, actual.fam);
    BOOST_CHECK_EQUAL(expected.child, actual.child);
    BOOST_CHECK_EQUAL(expected.dad, actual.dad);
    BOOST_CHECK_EQUAL(expected.mom, actual.mom);
    BOOST_CHECK(expected.sex == actual.sex);
    BOOST_CHECK_EQUAL(expected.sample_tree, actual.sample_tree);
}


BOOST_AUTO_TEST_CASE(test_skip_id) {

    std::vector<Pedigree::Member> expected_member_table {
        Pedigree::Member{0,0,0,0,Pedigree::Sex::Unknown,""},
        Pedigree::Member{0,1,0,0,Pedigree::Sex::Male,"NA12891"},
        Pedigree::Member{0,2,0,0,Pedigree::Sex::Female,"NA12892"},
        Pedigree::Member{0,3,0,0,Pedigree::Sex::Male,"NA12901"},
        Pedigree::Member{0,4,0,0,Pedigree::Sex::Female,"NA12902"},
        Pedigree::Member{0,5,1,2,Pedigree::Sex::Female,"NA12878"},
        Pedigree::Member{0,6,3,4,Pedigree::Sex::Female,"NA12903"}
    };

    std::vector<std::string> expected_name {
        "", "1", "2", "11", "12", "3", "13"
    };

    Pedigree ped;
    std::string input =
        "1\t1\t0\t0\t1\tNA12891\n"
        "1\t2\t0\t0\t2\tNA12892\n"
        "1\t3\t1\t2\t2\tNA12878\n"
        "1\t11\t0\t0\t1\tNA12901\n"
        "1\t12\t0\t0\t2\tNA12902\n"
        "1\t13\t11\t12\t2\tNA12903\n";

    ped.Parse(input);
    auto table = ped.table();
    BOOST_CHECK_EQUAL(expected_member_table.size(), table.size());
    BOOST_CHECK_EQUAL(expected_name.size(), ped.member_count());
    for (int i = 0; i < expected_name.size(); ++i) {
        boost_check_equal_member(expected_member_table[i], table[i]);
        auto name = ped.name(i);
        BOOST_CHECK_EQUAL(expected_name[i], name);
    }
}

BOOST_AUTO_TEST_CASE(test_multi_family_dup_id) {

    std::vector<Pedigree::Member> expected_member_table {
        Pedigree::Member{0,0,0,0,Pedigree::Sex::Unknown,""},
        Pedigree::Member{0,1,0,0,Pedigree::Sex::Male,"NA12891"},
        Pedigree::Member{0,2,0,0,Pedigree::Sex::Female,"NA12892"},
        Pedigree::Member{0,3,0,0,Pedigree::Sex::Male,"NA12901"},
        Pedigree::Member{0,4,0,0,Pedigree::Sex::Female,"NA12902"},
        Pedigree::Member{0,5,1,2,Pedigree::Sex::Female,"NA12878"},
        Pedigree::Member{0,6,3,4,Pedigree::Sex::Female,"NA12903"}
    };

    std::vector<std::string> expected_name {
        "", "1/1", "1/2", "2/1", "2/2", "1/3", "2/3"
    };

    Pedigree ped;
    std::string input =
        "1\t1\t0\t0\t1\tNA12891\n"
        "1\t2\t0\t0\t2\tNA12892\n"
        "1\t3\t1\t2\t2\tNA12878\n"
        "2\t1\t0\t0\t1\tNA12901\n"
        "2\t2\t0\t0\t2\tNA12902\n"
        "2\t3\t1\t2\t2\tNA12903\n";

    ped.Parse(input);
    auto table = ped.table();
    BOOST_CHECK_EQUAL(expected_member_table.size(), table.size());
    BOOST_CHECK_EQUAL(expected_name.size(), ped.member_count());
    for (int i = 0; i < expected_name.size(); ++i) {
        boost_check_equal_member(expected_member_table[i], table[i]);
        auto name = ped.name(i);
        BOOST_CHECK_EQUAL(expected_name[i], name);
    }
}



BOOST_AUTO_TEST_CASE(test_multi_family_diff_id) {

    std::vector<Pedigree::Member> expected_member_table {
        Pedigree::Member{0,0,0,0,Pedigree::Sex::Unknown,""},
        Pedigree::Member{0,1,0,0,Pedigree::Sex::Male,"NA12891"},
        Pedigree::Member{0,2,0,0,Pedigree::Sex::Female,"NA12892"},
        Pedigree::Member{0,3,0,0,Pedigree::Sex::Male,"NA12901"},
        Pedigree::Member{0,4,0,0,Pedigree::Sex::Female,"NA12902"},
        Pedigree::Member{0,5,1,2,Pedigree::Sex::Female,"NA12878"},
        Pedigree::Member{0,6,3,4,Pedigree::Sex::Female,"NA12903"}
    };

    std::vector<std::string> expected_name {
        "", "1", "2", "21", "22", "3", "23"
    };

    Pedigree ped;
    std::string input =
        "1\t1\t0\t0\t1\tNA12891\n"
        "1\t2\t0\t0\t2\tNA12892\n"
        "1\t3\t1\t2\t2\tNA12878\n"
        "2\t21\t0\t0\t1\tNA12901\n"
        "2\t22\t0\t0\t2\tNA12902\n"
        "2\t23\t21\t22\t2\tNA12903\n";

    ped.Parse(input);
    auto table = ped.table();
    BOOST_CHECK_EQUAL(expected_member_table.size(), table.size());
    BOOST_CHECK_EQUAL(expected_name.size(), ped.member_count());
    for (int i = 0; i < expected_name.size(); ++i) {
        boost_check_equal_member(expected_member_table[i], table[i]);
        auto name = ped.name(i);
        BOOST_CHECK_EQUAL(expected_name[i], name);
    }
}

BOOST_AUTO_TEST_CASE(test_dup_entries) {

    std::vector<Pedigree::Member> expected_member_table {
        Pedigree::Member{0,0,0,0,Pedigree::Sex::Unknown,""},
        Pedigree::Member{0,1,0,0,Pedigree::Sex::Male,"NA12891"},
        Pedigree::Member{0,2,0,0,Pedigree::Sex::Female,"NA12892"},
        Pedigree::Member{0,3,0,0,Pedigree::Sex::Male,"NA12901"},
        Pedigree::Member{0,4,0,0,Pedigree::Sex::Female,"NA12902"},
        Pedigree::Member{0,5,1,2,Pedigree::Sex::Female,"NA12878"},
        Pedigree::Member{0,6,3,4,Pedigree::Sex::Female,"NA12903"}
    };

    std::vector<std::string> expected_name {
        "", "1", "2", "21", "22", "3", "23"
    };

    Pedigree ped;
    std::string input =
        "1\t1\t0\t0\t1\tNA12891\n"
        "1\t2\t0\t0\t2\tNA12892\n"
        "1\t3\t1\t2\t2\tNA12878\n"
        "1\t3\t100\t200\t2\tNA12878\n" //dup
        "2\t21\t0\t0\t1\tNA12901\n"
        "2\t22\t0\t0\t2\tNA12902\n"
        "2\t22\t200\t200\t2\tNA12902\n" //dup
        "2\t23\t21\t22\t2\tNA12903\n";

    ped.Parse(input);
    auto table = ped.table();
    BOOST_CHECK_EQUAL(expected_member_table.size(), table.size());
    BOOST_CHECK_EQUAL(expected_name.size(), ped.member_count());
    for (int i = 0; i < expected_name.size(); ++i) {
        boost_check_equal_member(expected_member_table[i], table[i]);
        auto name = ped.name(i);
        BOOST_CHECK_EQUAL(expected_name[i], name);
    }
}


BOOST_AUTO_TEST_CASE(test_multi_family_dup_entries) {

    std::vector<Pedigree::Member> expected_member_table {
        Pedigree::Member{0,0,0,0,Pedigree::Sex::Unknown,""},
        Pedigree::Member{0,1,0,0,Pedigree::Sex::Male,"NA12891"},
        Pedigree::Member{0,2,0,0,Pedigree::Sex::Female,"NA12892"},
        Pedigree::Member{0,3,0,0,Pedigree::Sex::Male,"NA12901"},
        Pedigree::Member{0,4,0,0,Pedigree::Sex::Female,"NA12902"},
        Pedigree::Member{0,5,1,2,Pedigree::Sex::Female,"NA12878"},
        Pedigree::Member{0,6,3,4,Pedigree::Sex::Female,"NA12903"}
    };

    std::vector<std::string> expected_name {
        "", "1/1", "1/2", "2/1", "2/2", "1/3", "2/3"
    };

    Pedigree ped;
    std::string input =
        "1\t1\t0\t0\t1\tNA12891\n"
        "1\t2\t0\t0\t2\tNA12892\n"
        "1\t3\t1\t2\t2\tNA12878\n"
        "1\t3\t100\t200\t2\tNA12878\n" //dup
        "2\t1\t0\t0\t1\tNA12901\n"
        "2\t2\t0\t0\t2\tNA12902\n"
        "2\t2\t200\t200\t2\tNA12902\n" //dup
        "2\t3\t1\t2\t2\tNA12903\n";

    ped.Parse(input);
    auto table = ped.table();
    BOOST_CHECK_EQUAL(expected_member_table.size(), table.size());
    BOOST_CHECK_EQUAL(expected_name.size(), ped.member_count());
    for (int i = 0; i < expected_name.size(); ++i) {
        boost_check_equal_member(expected_member_table[i], table[i]);
        auto name = ped.name(i);
        BOOST_CHECK_EQUAL(expected_name[i], name);
    }
}


} // namespace io
} // namespace dng
