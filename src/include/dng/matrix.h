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

struct depth_t {
    int32_t counts[4]{0,0,0,0};
};

typedef std::vector<depth_t> RawDepths;

typedef Eigen::Array<double, Eigen::Dynamic, 1, 0, 10, 1> GenotypeArray;
typedef std::vector<GenotypeArray, Eigen::aligned_allocator<GenotypeArray>>
        GenotypeArrayVector;

#define DNG_INDIVIDUAL_BUFFER_MIN DBL_MIN
#define DNG_INDIVIDUAL_BUFFER_ONES GenotypeArrayVector::value_type::Ones(10)
#define DNG_INDIVIDUAL_BUFFER_ZEROS GenotypeArrayVector::value_type::Zero(10)

typedef Eigen::MatrixXd TransitionMatrix; // element (i,j) is the P(j|i)
typedef std::vector<TransitionMatrix> TransitionMatrixVector;

typedef Eigen::ArrayXXd TemporaryMatrix;

typedef Eigen::ArrayXd ParentArray;
typedef std::vector<ParentArray> ParentArrayVector;

typedef Eigen::Matrix4d MutationMatrix;

constexpr int folded_diploid_nucleotides[10][2] = {{0, 0}, {1, 1}, {2, 2}, {3, 3},
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};

// converts from unfolded genotype [0-15] to a folded genotype [0-10]
constexpr int folded_diploid_genotypes_matrix[4][4] = {
    {0, 4, 5, 6},
    {4, 1, 7, 8},
    {5, 7, 2, 9},
    {6, 8, 9, 3}};

constexpr int folded_diploid_genotypes[16] = {0, 4, 5, 6, 4, 1, 7, 8, 5, 7, 2, 9, 6, 8, 9, 3};

// converts from a folded genotype to an unfolded genotype
constexpr int unfolded_diploid_genotypes_upper[10] = {0, 5, 10, 15, 1, 2, 3, 6, 7, 11};
constexpr int unfolded_diploid_genotypes_lower[10] = {0, 5, 10, 15, 4, 8, 12, 9, 13, 14};

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
