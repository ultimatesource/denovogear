/*
 * Copyright (c) 2015-2016 Reed A. Cartwright
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
#include <utility>
#include <cstdint>
#include <cmath>
#include <cassert>
#include <initializer_list>

#include <dng/detail/unit_test.h>

namespace dng {
namespace stats {

double fisher_exact_test(int a11, int a12, int a21, int a22);

double g_test(double a11, double a12, double a21, double a22);

double ad_two_sample_test(std::vector<int> a, std::vector<int> b);

// Derived from Python's math.fsum
class ExactSum {
public:
    explicit ExactSum(double value=0.0) {
        partials_.reserve(32);
        if(value != 0.0) {
            partials_.push_back(value);
        }
    }
    ExactSum& operator=(double value) {
        partials_.clear();
        if(value != 0.0) {
            partials_.push_back(value);
        }
        return *this;
    }
    ExactSum& operator=(const ExactSum&) = default;
    ExactSum& operator=(ExactSum&&) = default;

    ExactSum& operator()(double x);
    ExactSum& operator()(const ExactSum& x);
    double result() const;

    ExactSum& operator+=(double value) {
        return operator()(value);
    }
    ExactSum& operator+=(const ExactSum& value) {
        return operator()(value);
    }

    ExactSum& add(double value) {
        return operator()(value);
    }
    ExactSum& add(const ExactSum& value) {
        return operator()(value);
    }

    operator double() const {
        return result();
    }
    bool operator==(double x) const {
        return (result() == x);
    }

    bool failed() const {
        return failed_;
    }
    operator bool() const {
        return failed_;
    }
    bool operator==(bool b) const {
        return failed_ == b;
    }

private:
    std::vector<double> partials_;
    double special_sum_{0.0};
    bool failed_{false};
};

template<typename SinglePassRange>
double exact_sum(SinglePassRange const & range) {
    ExactSum sum;
    for(auto &&x : range) {
        sum(x);
    }
    return sum;
}

template<typename T>
double exact_sum(const std::initializer_list<T>& range) {
    ExactSum sum;
    for(auto &&x : range) {
        sum(x);
    }
    return sum;
}

} // namespace stats
} // namespace dng

#endif // DNG_STATS_H
