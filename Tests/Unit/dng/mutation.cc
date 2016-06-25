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


#define BOOST_TEST_MODULE dng::mutation

#include <iostream>

#include <dng/mutation.h>

#include "../boost_test_helper.h"

using namespace dng;

struct CreateMutationMatrix {

    std::string fixture;

    MutationMatrix equal_mutation_matrix;
    MutationMatrix unequal_mutation_matrix;
    std::array<double, 4> equal_freq = {{0.25, 0.25, 0.25, 0.25}};
    std::array<double, 4> unequal_freq = {{0.1, 0.2, 0.3, 0.4}};


    CreateMutationMatrix(std::string s = "CreateMutationMatrix") : fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        equal_mutation_matrix = f81::matrix(1e-6, equal_freq);
        unequal_mutation_matrix = f81::matrix(1e-6, unequal_freq);

    }

    ~CreateMutationMatrix() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }


};

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }


BOOST_FIXTURE_TEST_SUITE(test_f81_suite, CreateMutationMatrix )

BOOST_AUTO_TEST_CASE(test_f81_equal) {

    auto diag = equal_mutation_matrix.diagonal();
    std::vector<double> expected_diag = {0.999999, 0.999999, 0.999999, 0.999999};
    std::vector<double> expected_off_diag = {3.333331111e-07, 3.333331111e-07, 3.333331111e-07, 3.333331111e-07};

    boost_check_close_vector(expected_diag, diag);
    for (int i = 0; i < 4; ++i) {
        auto actual = equal_mutation_matrix.col(i);
        Eigen::Array4d expected = Eigen::Array4d::Constant(expected_off_diag[i]);
        expected(i) = expected_diag[i];
        boost_check_close_vector(expected, actual);
    }

}

BOOST_AUTO_TEST_CASE(test_f81) {

    auto diag = unequal_mutation_matrix.diagonal();
    std::vector<double> expected_diag = {0.9999987143, 0.9999988571, 0.9999990000, 0.9999991429};
    std::vector<double> expected_off_diag = {1.428570408e-07, 2.857140816e-07, 4.285711225e-07, 5.714281633e-07};


    boost_check_close_vector(expected_diag, diag);
    for (int i = 0; i < 4; ++i) {
        auto actual = unequal_mutation_matrix.col(i);
        Eigen::Array4d expected = Eigen::Array4d::Constant(expected_off_diag[i]);
        expected(i) = expected_diag[i];
        boost_check_close_vector(expected, actual);
    }

}


BOOST_AUTO_TEST_CASE(test_f81_random) {

    auto matrix = f81::matrix(1e-3, {0.01, 0.1, 0.19, 0.7});
    auto diag = matrix.diagonal();
    std::vector<double> expected_diag = {0.9978677587, 0.9980615989, 0.9982554390, 0.9993538663};
    std::vector<double> expected_off_diag = {0.0000215377905, 0.0002153779050, 0.0004092180195, 0.0015076453352};


    boost_check_close_vector(expected_diag, diag);
    for (int i = 0; i < 4; ++i) {
        auto actual = matrix.col(i);
        Eigen::Array4d expected = Eigen::Array4d::Constant(expected_off_diag[i]);
        expected(i) = expected_diag[i];
        boost_check_close_vector(expected, actual);
    }

}



BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE(test_mutation_suite, CreateMutationMatrix )

BOOST_AUTO_TEST_CASE(test_mitosis_haploid_matrix) {


    TransitionMatrix expected_transition_matrix_negative = TransitionMatrix::Zero(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            expected_transition_matrix_negative(i, j) = unequal_mutation_matrix(i, j);
        }
    }
    TransitionMatrix expected_copy_unequal_mutation_matrix = unequal_mutation_matrix;


    TransitionMatrix expected_transition_matrix_zero = TransitionMatrix::Zero(4, 4);
    expected_transition_matrix_zero.diagonal() = unequal_mutation_matrix.diagonal();


    TransitionMatrix expected_transition_matrix_one = unequal_mutation_matrix;
    expected_transition_matrix_one.diagonal() = Eigen::Vector4d::Zero(4);


    auto actual_negative = mitosis_haploid_matrix(unequal_mutation_matrix, -1);
    boost_check_matrix(expected_transition_matrix_negative, actual_negative);
    boost_check_matrix(expected_copy_unequal_mutation_matrix, actual_negative);

    auto actual_zero = mitosis_haploid_matrix(unequal_mutation_matrix, 0);
    boost_check_matrix(expected_transition_matrix_zero, actual_zero);

    auto actual_one = mitosis_haploid_matrix(unequal_mutation_matrix, 1);
    boost_check_matrix(expected_transition_matrix_one, actual_one);

}

BOOST_AUTO_TEST_CASE(test_mitosis_diploid_matrix) {

    TransitionMatrix ret = TransitionMatrix::Zero(10, 10);
//    for(int i = 0; i < 10; ++i) {
//        int h = unfolded_diploid_genotypes_upper[i];
//        for(int j = 0; j < 16; ++j) {
//            int k = folded_diploid_genotypes[j];
//            ret(i, k) += (mutype < 0 || mutype == mitotic_diploid_mutation_counts[h][j]) ?
//                         kronecker_product(m, m, h, j) : 0.0;
//        }
//    }
    auto kp = kroneckerProduct(unequal_mutation_matrix, unequal_mutation_matrix);

    TransitionMatrix expected_transition_matrix_negative = TransitionMatrix::Zero(10, 10);
    for(int m = 0; m < 10; ++m) {
        int i = unfolded_diploid_genotypes_upper[m];

    for (int j = 0; j < 16; ++j) {
            int n = folded_diploid_genotypes[j];
            expected_transition_matrix_negative(m, n) += kp.coeff(i, j);
        }
    }

    TransitionMatrix expected_transition_matrix_zero = TransitionMatrix::Zero(10, 10);

    TransitionMatrix expected_transition_matrix_one = TransitionMatrix::Zero(10, 10);

    TransitionMatrix expected_transition_matrix_two = TransitionMatrix::Zero(10, 10);


    auto actual_negative = mitosis_diploid_matrix(unequal_mutation_matrix, -1);
    boost_check_matrix(expected_transition_matrix_negative, actual_negative);

    auto actual_zero = mitosis_diploid_matrix(unequal_mutation_matrix, 0);
//    boost_check_matrix(expected_transition_matrix_zero, actual_zero);

    auto actual_one = mitosis_diploid_matrix(unequal_mutation_matrix, 1);
//    boost_check_matrix(expected_transition_matrix_one, actual_one);

    auto actual_two = mitosis_diploid_matrix(unequal_mutation_matrix, 2);
//    boost_check_matrix(expected_transition_matrix_two, actual_two);


}

//TODO(SHW): Implement these tests
BOOST_AUTO_TEST_CASE(test_meiosis_haploid_matrix) {
}
BOOST_AUTO_TEST_CASE(test_meiosis_diploid_matrix) {
}
BOOST_AUTO_TEST_CASE(test_mitosis_haploid_mean_matrix) {
}
BOOST_AUTO_TEST_CASE(test_mitosis_diploid_mean_matrix) {
}
BOOST_AUTO_TEST_CASE(test_meiosis_haploid_mean_matrix) {
}
BOOST_AUTO_TEST_CASE(test_meiosis_diploid_mean_matrix) {
}

BOOST_AUTO_TEST_SUITE_END()

