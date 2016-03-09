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

#define BOOST_TEST_MODULE dng::lib::peeling

#include <boost/test/unit_test.hpp>
//#include <boost/test/tools/floating_point_comparison.hpp>

#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <dng/peeling.h>

#include <boost_test_helper.h>
#include <fixture/fixture_random_family.h>

//#include <boost/test/data/test_case.hpp>
////#include <boost/test/data/monomorphic.hpp>
//namespace data = boost::unit_test::data;
// TODO: Example of BOOST_DATA_TEST_CASE and BOOST_PARAM_TEST_CASE.
// TODO: Should be able to replace the for loop with these.
// TODO: Might not be able to use fixture.
//
//
//std::vector<int> test_types;//
//
//BOOST_AUTO_TEST_SUITE(suite1,
//  * utf::fixture<RandomFamily>(std::string("FX"))
//  * utf::fixture<RandomFamily>(std::string("FX2")))
//
//  BOOST_AUTO_TEST_CASE(test1, * utf::fixture(&setup, &teardown))
//  {
//    BOOST_TEST_MESSAGE("running test1");
//    BOOST_TEST(true);
//  }
//
//  BOOST_AUTO_TEST_CASE(test2)
//  {
//    BOOST_TEST_MESSAGE("running test2");
//    BOOST_TEST(true);
//  }
//
//BOOST_DATA_TEST_CASE( test_case_arity1, data::xrange(5), my_var )
//{
//    BOOST_TEST_MESSAGE("running data: ");
//    BOOST_TEST((my_var <= 4 && my_var >= 0));
//}
//BOOST_PARAM_TEST_CASE(test_function, params_begin, params_end);
//BOOST_AUTO_TEST_SUITE_END()


using namespace dng;
namespace utf = boost::unit_test;

const int NUM_TEST = 100;


void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }



//BOOST_AUTO_TEST_SUITE(test_peeling_suite,  * utf::fixture<RandomFamily>(std::string("FX")) )
BOOST_FIXTURE_TEST_SUITE(test_peeling_suite, RandomFamily)

    BOOST_AUTO_TEST_CASE(test_sum_over_child, *utf::fixture(&setup, &teardown)) {
        for (int t = 0; t <NUM_TEST; ++t) {

            init_family();

            PairedGenotypeArray expected = PairedGenotypeArray::Ones(100, 1);
            for (int k = 2; k < total_family_size; ++k) {
                PairedGenotypeArray temp_array = PairedGenotypeArray::Zero(100, 1);

                for (int i = 0; i < trans_matrix[k].rows(); ++i) {
                    double x = 0;
                    for (int j = 0; j < trans_matrix[k].cols(); ++j) {
                        x += trans_matrix[k](i, j) * lower_array[k][j];
                    }
                    temp_array(i, 0) = x;
                }
                expected *= temp_array;
            }
            //Done expected


            dng::peel::workspace_t workspace;
            dng::TransitionVector full_matrix;
            copy_family_to_workspace(workspace, full_matrix, total_family_size,
                                     lower_array, upper_array,
                                     trans_matrix);

            PairedGenotypeArray result = dng::peel::sum_over_children(workspace, family, full_matrix);

            BoostCheckMatrix(expected, result, 100, 1);
        }
    }


    BOOST_AUTO_TEST_CASE(test_up_core, *utf::fixture(&setup, &teardown)) {

        for (int t = 0; t <NUM_TEST; ++t) {

            init_family_parent_child_only();

            GenotypeArray expected = GenotypeArray::Zero();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    expected[i] += trans_matrix[1](i, j) * lower_array[1][j];
                }
            }
            //Done expected

            dng::peel::workspace_t workspace;
            dng::TransitionVector full_matrix;
            copy_family_to_workspace(workspace, full_matrix, total_family_size,
                                     lower_array, upper_array, trans_matrix);

            GenotypeArray result = dng::peel::up_core(workspace, family, full_matrix);
            BoostCheckCloseVector(expected, result, 10);

        }
    }



    BOOST_AUTO_TEST_CASE(test_down) {

        for (int t = 0; t < NUM_TEST; ++t) {
            init_family_parent_child_only();

            GenotypeArray expected = GenotypeArray::Zero();
            GenotypeArray expected_fast = GenotypeArray::Zero();

            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    expected_fast[i] += trans_matrix[1](j, i) * upper_array[0][j];
                    expected[i] += trans_matrix[1](j, i) * upper_array[0][j] * lower_array[0][j];
                }
            }
            //Done expected

            dng::peel::workspace_t workspace;
            dng::TransitionVector full_matrix;
            copy_family_to_workspace(workspace, full_matrix, total_family_size,
                                     lower_array, upper_array, trans_matrix);

            dng::peel::down(workspace, family, full_matrix);
            GenotypeArray result = workspace.upper[1];

            dng::peel::down_fast(workspace, family, full_matrix);
            GenotypeArray result_fast = workspace.upper[1];

            BoostCheckCloseVector(expected, result, 10);
            BoostCheckCloseVector(expected_fast, result_fast, 10);

        }
    }

    BOOST_AUTO_TEST_CASE(test_up) {

        for (int t = 0; t < NUM_TEST; ++t) {
            init_family_parent_child_only();

            GenotypeArray expected = GenotypeArray::Zero();
            GenotypeArray expected_fast = GenotypeArray::Zero();

            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    expected[i] += trans_matrix[1](i, j) * lower_array[1][j];
                }
                expected_fast[i] = expected[i];
                expected[i] *= lower_array[0][i];
            }
            //Done expected

            dng::peel::workspace_t workspace;
            dng::TransitionVector full_matrix;
            copy_family_to_workspace(workspace, full_matrix, total_family_size,
                                     lower_array, upper_array, trans_matrix);

            dng::peel::up(workspace, family, full_matrix);
            GenotypeArray result = workspace.lower[0];

            dng::peel::up_fast(workspace, family, full_matrix);
            GenotypeArray result_fast = workspace.lower[0];

            BoostCheckCloseVector(expected, result, 10);
            BoostCheckCloseVector(expected_fast, result_fast, 10);

        }
    }


    BOOST_AUTO_TEST_CASE(test_to_father) {

        for (int t = 0; t < NUM_TEST; ++t) {

            init_family();

            PairedGenotypeArray all_child = PairedGenotypeArray::Ones(100, 1);
            for (int c = CHILD_OFFSET; c < total_family_size; ++c) {
                PairedGenotypeArray one_child = PairedGenotypeArray::Zero(100, 1);
                for (int i = 0; i < 100; ++i) {
                    for (int j = 0; j < 10; ++j) {
                        one_child(i, 0) += trans_matrix[c](i, j) * lower_array[c][j];
                    }
                }
                all_child *= one_child;
            }
            all_child.resize(10, 10);

            GenotypeArray ga_mum;
            for (int j = 0; j < 10; ++j) {
                ga_mum[j] = upper_array[1][j] * lower_array[1][j];
            }

            GenotypeArray expected = GenotypeArray::Zero();
            GenotypeArray expected_fast = GenotypeArray::Zero();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    expected[i] += all_child(i, j) * ga_mum[j];
                }
                expected_fast[i] = expected[i];
                expected[i] *= lower_array[0][i];
            }
            //Done expected

            dng::peel::workspace_t workspace;
            dng::TransitionVector full_matrix;
            copy_family_to_workspace(workspace, full_matrix, total_family_size,
                                     lower_array, upper_array, trans_matrix);


            dng::peel::to_father(workspace, family, full_matrix);
            GenotypeArray result = workspace.lower[0];
            for (int i = 0; i < 10; ++i) {
                BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
            }

            dng::peel::to_father_fast(workspace, family, full_matrix);
            GenotypeArray result_fast = workspace.lower[0];


            BoostCheckCloseVector(expected, result, 10);
            BoostCheckCloseVector(expected_fast, result_fast, 10);



        }
    }


    BOOST_AUTO_TEST_CASE(test_to_mother) {

        for (int t = 0; t < NUM_TEST; ++t) {

            init_family();

            PairedGenotypeArray all_child = PairedGenotypeArray::Ones(100, 1);
            for (int c = CHILD_OFFSET; c < total_family_size; ++c) {
                PairedGenotypeArray one_child = PairedGenotypeArray::Zero(100, 1);
                for (int i = 0; i < 100; ++i) {
                    for (int j = 0; j < 10; ++j) {
                        one_child(i, 0) += trans_matrix[c](i, j) * lower_array[c][j];
                    }
                }
                all_child *= one_child;
            }
            all_child.resize(10, 10);

            GenotypeArray ga_dad;
            for (int j = 0; j < 10; ++j) {
                ga_dad[j] = upper_array[0][j] * lower_array[0][j];
            }

            GenotypeArray expected = GenotypeArray::Zero();
            GenotypeArray expected_fast = GenotypeArray::Zero();
            auto temp_child = all_child.transpose();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    expected[i] += temp_child(i, j) * ga_dad[j];
                }
                expected_fast[i] = expected[i];
                expected[i] *= lower_array[1][i];
            }
            //Done expected

            dng::peel::workspace_t workspace;
            dng::TransitionVector full_matrix;
            copy_family_to_workspace(workspace, full_matrix, total_family_size,
                                     lower_array, upper_array, trans_matrix);


            dng::peel::to_mother(workspace, family, full_matrix);
            GenotypeArray result = workspace.lower[1];

            dng::peel::to_mother_fast(workspace, family, full_matrix);
            GenotypeArray result_fast = workspace.lower[1];

            BoostCheckCloseVector(expected, result, 10);
            BoostCheckCloseVector(expected_fast, result_fast, 10);


        }
    }


    BOOST_AUTO_TEST_CASE(test_to_child) {

        rand_unif = std::uniform_int_distribution<>(2,10); //min 2 childred required

        for (int t = 0; t < NUM_TEST; ++t) {

            init_family();

            PairedGenotypeArray all_child = PairedGenotypeArray::Ones(100, 1);
            for (int c = CHILD_OFFSET + 1; c < total_family_size; ++c) {
                PairedGenotypeArray one_child = PairedGenotypeArray::Zero(100, 1);
                for (int i = 0; i < 100; ++i) {
                    for (int j = 0; j < 10; ++j) {
                        one_child(i, 0) += trans_matrix[c](i, j) * lower_array[c][j];
                    }
                }
                all_child *= one_child;
            }

            GenotypeArray dad_array;
            GenotypeArray mom_array;
            for (int j = 0; j < 10; ++j) {
                dad_array[j] = upper_array[0][j] * lower_array[0][j];
                mom_array[j] = upper_array[1][j] * lower_array[1][j];
            }
            PairedGenotypeArray parents_array = PairedGenotypeArray::Zero(100, 1);
            for (int d = 0; d < 10; ++d) {
                for (int m = 0; m < 10; ++m) {
                    parents_array(d * 10 + m, 0) = dad_array[d] * mom_array[m];
                }
            }
            PairedGenotypeArray parents_children = parents_array * all_child;

            GenotypeArray expected = GenotypeArray::Zero();
            GenotypeArray expected_fast = GenotypeArray::Zero();
            auto temp_m = trans_matrix[2].transpose();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 100; ++j) {
                    expected[i] += temp_m(i, j) * parents_children(j, 0);
                    expected_fast[i] += temp_m(i, j) * parents_array(j, 0);
                }
            }
            //Done expected

            dng::peel::workspace_t workspace;
            dng::TransitionVector full_matrix;
            copy_family_to_workspace(workspace, full_matrix, total_family_size,
                                     lower_array, upper_array, trans_matrix);


            dng::peel::to_child(workspace, family, full_matrix);
            GenotypeArray result = workspace.upper[2];

            peel::family_members_t family_fast{family[0], family[1], family[2]};
            dng::peel::to_child_fast(workspace, family_fast, full_matrix);
            GenotypeArray result_fast = workspace.upper[2];


            BoostCheckCloseVector(expected, result, 10);
            BoostCheckCloseVector(expected_fast, result_fast, 10);

        }
    }

BOOST_AUTO_TEST_SUITE_END()

