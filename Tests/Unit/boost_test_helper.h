//
// Created by steven on 1/15/16.
//
// Replace some with with BOOST later
//

#ifndef DENOVOGEAR_BOOST_TEST_HELPER_H
#define DENOVOGEAR_BOOST_TEST_HELPER_H

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

//TODO: Where should this file go?
//FIXME: too many global
//TODO(Reed): I think we should make this configurable somehow from a CMake run. It is important to test both absolute equality and threshold-based equality.

const double ABS_TEST_THRESHOLD = 1e-10;
const double BOOST_CLOSE_THRESHOLD = 1e-6;


template<typename A>
void boost_check_array(A &expected, A &result, int expected_size) {
    
    BOOST_CHECK_EQUAL(expected_size, expected.size());
    BOOST_CHECK_EQUAL(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
    }

}

template<typename M>
void boost_check_matrix(M &expected, M &result) {
	BOOST_CHECK_EQUAL(expected.rows(), result.rows());
	BOOST_CHECK_EQUAL(expected.cols(), result.cols());

	for (int i = 0; i < expected.rows(); i++) {
		for (int j = 0; j < expected.cols(); j++) {
			if (expected(i, j) == 0) {
				BOOST_CHECK_EQUAL(0, result(i,j));
			} else {
				BOOST_CHECK_CLOSE(expected(i, j), result(i, j),
						BOOST_CLOSE_THRESHOLD);
			}
		}
	}
}


template<typename M>
void boost_check_matrix(M &expected, M &result, int expected_rows,
		int expected_cols) {
	BOOST_CHECK_EQUAL(expected_rows, expected.rows());
	BOOST_CHECK_EQUAL(expected_cols, expected.cols());

	boost_check_matrix(expected, result);
}

template<typename V, typename V2>
void boost_check_close_vector(V &expected, V2 &result) {

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

template<typename A, typename B>
void AssertTrue(A expected, B actual){
    assert(expected==actual);
};

template<typename A>
void AssertNear(A expected, A actual){
    assert(((expected - actual)/expected) < ABS_TEST_THRESHOLD);

};

template<typename A>
void AssertEigenMatrixNear(A expected, A actual){
    for (int j = 0; j < expected.rows(); ++j) {
        for (int k = 0; k < expected.cols(); ++k) {
            AssertNear(expected(j,k), actual(j,k));
        }
    }
};


#endif //DENOVOGEAR_BOOST_TEST_HELPER_H
