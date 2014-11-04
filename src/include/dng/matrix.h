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

	typedef Eigen::MatrixXd TransitionMatrix;
	typedef std::vector<TransitionMatrix> TransitionVector;

	typedef Eigen::ArrayXd GenotypeArray;
	typedef Eigen::ArrayXXd PairedGenotypeArray;
	typedef std::vector<GenotypeArray> IndividualBuffer;


	//typedef std::vector<GenotypeArray, Eigen::aligned_allocator<GenotypeArray>> IndividualBuffer;

	const int nucleotides[10][2] = {{0,0},{0,1},{0,2},{0,3},
		{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}};
};

#endif // DNG_MATRIX_H
