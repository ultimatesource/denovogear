//
// Created by steven on 1/29/16.
//

#ifndef DENOVOGEAR_ASSERT_UTILS_H
#define DENOVOGEAR_ASSERT_UTILS_H

#include <cassert>


const double ASSERT_CLOSE_THRESHOLD = 1e-6;

//TODO: HACK: FIXME:
template<typename A, typename B>
void AssertEqual(A expected, B actual) {
    assert(expected == actual);
};

template<typename A>
void AssertEqual(A expected, A actual) {
    assert(expected == actual);
};

template<typename A>
void AssertNear(A expected, A actual) {
    if (!(expected == 0 && actual == 0)) {
        assert(((expected - actual) / expected) < ASSERT_CLOSE_THRESHOLD);
    }

};

//template<typename A>
//void AssertVectorNear(A expected, A actual){
//    AssertEqual(expected.size(), actual.size());
//    for (int j = 0; j < expected.size(); ++j) {
//        AssertNear(expected[j], actual[j]);
//    }
//};

template<typename A, typename B>
void AssertVectorEqual(A expected, B actual) {
    AssertEqual(expected.size(), actual.size());
    for (int j = 0; j < expected.size(); ++j) {
        AssertEqual(expected[j], actual[j]);
    }
};


template<typename A, typename B>
void AssertVectorNear(A expected, B actual) {
    AssertEqual(expected.size(), actual.size());
    for (int j = 0; j < expected.size(); ++j) {
        AssertNear(expected[j], actual[j]);
    }
};


template<typename A>
void AssertEigenMatrixNear(A expected, A actual) {
    AssertEqual(expected.rows(), actual.rows());
    AssertEqual(expected.cols(), actual.cols());
    for (int j = 0; j < expected.rows(); ++j) {
        for (int k = 0; k < expected.cols(); ++k) {
            AssertNear(expected(j, k), actual(j, k));
        }
    }
}


#endif //DENOVOGEAR_ASSERT_UTILS_H
