/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2017 Reed A. Cartwright
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
 *           Reed A. Cartwright <reed@cartwrig.ht>
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

#define BOOST_TEST_MODULE dng::relationship_graph

#include <dng/relationship_graph.h>

#include "../testing.h"

using namespace dng;
using namespace dng::detail;
using Sex = dng::Pedigree::Sex;

struct topology_t {
    std::string name;
    libraries_t libs;
    Pedigree ped;
};

struct rates_t {
    double germline;
    double somatic;
    double library;
};

template<typename F>
void run_graph_tests(F test, double prec = 2*DBL_EPSILON) {
    topology_t trio{"trio"};
    trio.libs = {
        {"Mom", "Dad", "Eve"},
        {"Mom", "Dad", "Eve"}
    };
    trio.ped.AddMember("Dad","0","0",Sex::Male,"");
    trio.ped.AddMember("Mom","0","0",Sex::Female,"");
    trio.ped.AddMember("Eve","Dad","Mom",Sex::Female,"");

    test(trio, prec);
}

BOOST_AUTO_TEST_CASE(test_inheritance_model_from_string) {
    BOOST_CHECK(inheritance_model("autosomal") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("maternal") == InheritanceModel::Maternal);
    BOOST_CHECK(inheritance_model("paternal") == InheritanceModel::Paternal);
    BOOST_CHECK(inheritance_model("x-linked") == InheritanceModel::XLinked);
    BOOST_CHECK(inheritance_model("y-linked") == InheritanceModel::YLinked);
    BOOST_CHECK(inheritance_model("w-linked") == InheritanceModel::WLinked);
    BOOST_CHECK(inheritance_model("z-linked") == InheritanceModel::ZLinked);
    
    BOOST_CHECK(inheritance_model("mitochondrial") == InheritanceModel::Maternal);
    BOOST_CHECK(inheritance_model("xlinked") == InheritanceModel::XLinked);
    BOOST_CHECK(inheritance_model("ylinked") == InheritanceModel::YLinked);
    BOOST_CHECK(inheritance_model("wlinked") == InheritanceModel::WLinked);
    BOOST_CHECK(inheritance_model("zlinked") == InheritanceModel::ZLinked);

    BOOST_CHECK(inheritance_model("a") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("m") == InheritanceModel::Maternal);
    BOOST_CHECK(inheritance_model("p") == InheritanceModel::Paternal);
    BOOST_CHECK(inheritance_model("x") == InheritanceModel::XLinked);
    BOOST_CHECK(inheritance_model("y") == InheritanceModel::YLinked);
    BOOST_CHECK(inheritance_model("w") == InheritanceModel::WLinked);
    BOOST_CHECK(inheritance_model("z") == InheritanceModel::ZLinked);

    BOOST_CHECK(inheritance_model("A") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("M") == InheritanceModel::Maternal);
    BOOST_CHECK(inheritance_model("P") == InheritanceModel::Paternal);
    BOOST_CHECK(inheritance_model("X") == InheritanceModel::XLinked);
    BOOST_CHECK(inheritance_model("Y") == InheritanceModel::YLinked);
    BOOST_CHECK(inheritance_model("W") == InheritanceModel::WLinked);
    BOOST_CHECK(inheritance_model("Z") == InheritanceModel::ZLinked);

    BOOST_CHECK(inheritance_model("Autosomal") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("AUTOSOMAL") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("au") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("aut") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("auto") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("autos") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("autoso") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("autosom") == InheritanceModel::Autosomal);
    BOOST_CHECK(inheritance_model("autosoma") == InheritanceModel::Autosomal);

    BOOST_CHECK_THROW(inheritance_model(""), std::invalid_argument);
    BOOST_CHECK_THROW(inheritance_model("unknown"), std::invalid_argument);
    BOOST_CHECK_THROW(inheritance_model("aut0somal"), std::invalid_argument);
    BOOST_CHECK_THROW(inheritance_model("b"), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_inheritance_model_to_string) {
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::Autosomal), "AUTOSOMAL");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::Maternal), "MATERNAL");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::Paternal), "PATERNAL");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::XLinked), "X-LINKED");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::YLinked), "Y-LINKED");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::WLinked), "W-LINKED");
    BOOST_CHECK_EQUAL(to_string(InheritanceModel::ZLinked), "Z-LINKED");

    BOOST_CHECK_THROW(to_string(InheritanceModel::Unknown), std::invalid_argument);
    BOOST_CHECK_THROW(to_string(static_cast<InheritanceModel>(-1)), std::invalid_argument);
}
