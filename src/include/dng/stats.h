/*
 * Copyright (c) 2015 Reed A. Cartwright
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
#ifndef DNG_STATS_H
#define DNG_STATS_H

#include <vector>
#include <cstdint>

namespace dng {
namespace stats {

double fisher_exact_test(int a11, int a12, int a21, int a22);

double g_test(double a11, double a12, double a21, double a22);

double ad_two_sample_test(std::vector<uint8_t> a, std::vector<uint8_t> b);

} // namespace stats
} // namespace dng

#endif // DNG_STATS_H
