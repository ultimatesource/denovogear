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

	typedef Eigen::Array<double, 10, 1> Vector10d;
	typedef Eigen::Matrix<double, 10, 10> Matrix10d;
	typedef Eigen::Matrix<double, 10, 10, Eigen::RowMajor> RowMatrix10d;
	typedef Eigen::Array<double, 100, 1> Vector100d;
	typedef Eigen::Matrix<double, 100, 100> Matrix100d;
	
	typedef Eigen::Matrix<double, 100, 10> MeiosisMatrix;
	typedef Matrix10d MitosisMatrix;

	typedef std::vector<Vector10d, Eigen::aligned_allocator<Vector10d>> IndividualBuffer;
	typedef std::vector<MeiosisMatrix, Eigen::aligned_allocator<MeiosisMatrix>> MeiosisMatrixVector;
	typedef std::vector<Matrix10d, Eigen::aligned_allocator<Matrix10d>> MitosisMatrixVector;

	const int nucleotides[10][2] = {{0,0},{0,1},{0,2},{0,3},
		{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}};
};

#endif // DNG_MATRIX_H
