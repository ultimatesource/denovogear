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


#ifndef DENOVOGEAR_BOOST_TEST_HELPER_H
#define DENOVOGEAR_BOOST_TEST_HELPER_H

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

//TODO: Where should this file go?
//FIXME: too many global
const double ABS_TEST_THRESHOLD = 1e-10;
const double BOOST_CLOSE_THRESHOLD = 1e-3;//
const double BOOST_SMALL_THRESHOLD = sqrt(std::numeric_limits<double>::epsilon());
//TODO: What is a good cut? 1e-8?// sqrt(std::numeric_limits<double>::epsilon());


template<typename V, typename V2>
void BoostCheckEqualVector(V &expected, V2 &result) {

    BOOST_CHECK_EQUAL(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        BOOST_CHECK(expected[i] == result[i]);
    }
}

template<typename V, typename V2>
void BoostCheckCloseVector(V &expected, V2 &result) {

    BOOST_CHECK_EQUAL(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        if(expected[i] == 0){
            BOOST_CHECK_EQUAL(0, result[i]);
        }
        else {
            BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
        }
    }
}

template<typename V, typename V2>
void BoostCheckCloseVector(V &&expected, V2 &result) {

    BOOST_CHECK_EQUAL(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        if(expected[i] == 0){
            BOOST_CHECK_EQUAL(0, result[i]);
        }
        else {
            BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
        }
    }
}

template<typename V, typename V2>
void BoostCheckCloseVector(V &expected, V2 &result, int expected_size) {
    
    BOOST_CHECK_EQUAL(expected_size, expected.size());
    BoostCheckCloseVector(expected, result);
}



template<typename M>
void BoostCheckMatrix(M &expected, M &result) {
    BOOST_CHECK_EQUAL(expected.rows(), result.rows());
    BOOST_CHECK_EQUAL(expected.cols(), result.cols());

    for (int i = 0; i < expected.rows(); ++i) {
        for (int j = 0; j < expected.cols(); ++j) {
            if(expected(i,j) == 0 ){
                BOOST_CHECK_EQUAL(0, result(i,j));
            }
            else {
                BOOST_CHECK_CLOSE(expected(i, j), result(i, j), BOOST_CLOSE_THRESHOLD);
            }
        }
    }
}


template<typename M>
void BoostCheckMatrix(M &expected, M &result, int expected_rows, int expected_cols) {
    BOOST_CHECK_EQUAL(expected_rows, expected.rows());
    BOOST_CHECK_EQUAL(expected_cols, expected.cols());

    BoostCheckMatrix(expected, result);
}




#endif //DENOVOGEAR_BOOST_TEST_HELPER_H
