/*
 * Copyright (c) 2014-2016 Reed A. Cartwright
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
#pragma once
#ifndef DNG_MATRIX_H
#define DNG_MATRIX_H

#include <cstdint>
#include <cfloat>
#include <algorithm>
#include <initializer_list>

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

#include <dng/detail/unit_test.h>

namespace dng {

typedef Eigen::ArrayXd GenotypeArray;
typedef std::vector<GenotypeArray, Eigen::aligned_allocator<GenotypeArray>>
        GenotypeArrayVector;

#define DNG_INDIVIDUAL_BUFFER_MIN DBL_MIN

typedef Eigen::MatrixXd TransitionMatrix; // element (i,j) is the P(j|i)
typedef std::vector<TransitionMatrix> TransitionMatrixVector;

typedef Eigen::ArrayXXd TemporaryMatrix;

typedef Eigen::ArrayXd ParentArray;
typedef std::vector<ParentArray> ParentArrayVector;

template<typename A, typename B>
inline auto kronecker_product_coef(const A &a, const B &b, std::size_t i,
                              std::size_t j) -> decltype(a(0, 0)*b(0, 0)) {
    assert(i < a.rows()*b.rows() && j < a.cols()*b.cols());
    return a(i / b.rows(), j / b.cols()) * b(i % b.rows(), j % b.cols());
}

// From http://stackoverflow.com/a/16287999
template <typename T, int R, int C>
inline T sum_kahan(const Eigen::Array<T, R, C> &xs) {
    if(xs.size() == 0) {
        return 0;
    }
    T sumP(0), sumN(0), tP(0), tN(0);
    T cP(0), cN(0), yP(0), yN(0);
    for(size_t i = 0, size = xs.size(); i < size; i++) {
        T temporary = (*(xs.data() + i));
        if(temporary > 0) {
            yP = temporary - cP;
            tP = sumP + yP;
            cP = (tP - sumP) - yP;
            sumP = tP;
        } else {
            yN = temporary - cN;
            tN = sumN + yN;
            cN = (tN - sumN) - yN;
            sumN = tN;
        }
    }
    return sumP + sumN;
}

} // namespace dng

#endif // DNG_MATRIX_H
