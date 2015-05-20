/*
 * Copyright (c) 2014 Reed A. Cartwright
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

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/KroneckerProduct>

namespace dng {

union depth_t {
    uint64_t key;
    uint16_t counts[4];
};

union depth5_t {
    struct {
        uint64_t key;
        uint64_t pad;
    };
    uint16_t counts[8];
};

#ifdef DNG_USE_DYNAMIC_GENOTYPE_ARRAY
typedef Eigen::ArrayXd GenotypeArray;
typedef std::vector<GenotypeArray> IndividualBuffer;
#else
typedef Eigen::Array<double, 10, 1> GenotypeArray;
typedef std::vector<GenotypeArray, Eigen::aligned_allocator<GenotypeArray>>
        IndividualBuffer;
#endif

#define DNG_INDIVIDUAL_BUFFER_ASSIGN_TYPE IndividualBuffer::value_type::Ones(10)

typedef Eigen::MatrixXd TransitionMatrix;
typedef std::vector<TransitionMatrix> TransitionVector;

typedef Eigen::ArrayXXd PairedGenotypeArray;

constexpr int nucleotides[10][2] = {{0, 0}, {0, 1}, {0, 2}, {0, 3},
    {1, 1}, {1, 2}, {1, 3}, {2, 2}, {2, 3}, {3, 3}
};

constexpr int genotypes_folded[16] = {0, 1, 2, 3, 1, 4, 5, 6, 2, 5, 7, 8, 3, 6, 8, 9};
};

#endif // DNG_MATRIX_H
