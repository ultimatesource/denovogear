//
// Created by steven on 1/15/16.
//
// Replace some with with BOOST later
//

#ifndef DENOVOGEAR_ASSERT_HELPER_H
#define DENOVOGEAR_ASSERT_HELPER_H

#include <boost/test/test_tools.hpp>

//FIXME: too many global
const double TEST_THRESHOLD = 1e-10;


template<typename A, typename B>
void AssertTrue(A expected, B actual){
    assert(expected==actual);
};

template<typename A>
void AssertNear(A expected, A actual){
    assert( ((expected - actual)/expected) < TEST_THRESHOLD);

};

template<typename A>
void AssertEigenMatrixNear(A expected, A actual){
    for (int j = 0; j < expected.rows(); ++j) {
        for (int k = 0; k < expected.cols(); ++k) {
            AssertNear(expected(j,k), actual(j,k));
        }
    }
};



#endif //DENOVOGEAR_ASSERT_HELPER_H
