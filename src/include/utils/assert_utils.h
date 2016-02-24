//
// Created by steven on 1/29/16.
//

#ifndef DENOVOGEAR_ASSERT_UTILS_H
#define DENOVOGEAR_ASSERT_UTILS_H

#include <cassert>


const double ASSERT_CLOSE_THRESHOLD = 1e-6;


//TODO: HACK: FIXME:
template<typename A, typename B>
void assert_equal(A expected, B actual) {
    assert(expected == actual);
};

template<typename A>
void assert_equal(A expected, A actual) {
    assert(expected == actual);
};

template<typename A>
void assert_near(A expected, A actual) {

#ifndef	NDEBUG
    if (!(expected == 0 && actual == 0)) {
        assert(((expected - actual) / expected) < ASSERT_CLOSE_THRESHOLD);
    }
#endif
};

//template<typename A>
//void AssertVectorNear(A expected, A actual){
//    assert_equal(expected.size(), actual.size());
//    for (int j = 0; j < expected.size(); ++j) {
//        AssertNear(expected[j], actual[j]);
//    }
//};

template<typename A, typename B>
void assert_equal_vector(A expected, B actual) {
#ifndef	NDEBUG
    assert_equal(expected.size(), actual.size());
    for (int j = 0; j < expected.size(); ++j) {
        assert_equal(expected[j], actual[j]);
    }
#endif
};


template<typename A, typename B>
void assert_near_vector(A expected, B actual) {
#ifndef	NDEBUG
    assert_equal(expected.size(), actual.size());
    for (int j = 0; j < expected.size(); ++j) {
        assert_near(expected[j], actual[j]);
    }
#endif
};


template<typename A>
void assert_near_eigen_matrix(A expected, A actual) {
#ifndef	NDEBUG
    assert_equal(expected.rows(), actual.rows());
    assert_equal(expected.cols(), actual.cols());
    for (int j = 0; j < expected.rows(); ++j) {
        for (int k = 0; k < expected.cols(); ++k) {
            assert_near(expected(j, k), actual(j, k));
        }
    }
#endif
}


#endif //DENOVOGEAR_ASSERT_UTILS_H
