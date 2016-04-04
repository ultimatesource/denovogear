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

double ad_two_sample_test(std::vector<uint8_t> a, std::vector<uint8_t> b);

// Derived from Python's math.fsum
class ExactSum {
public:
    inline ExactSum(double value=0.0) {
        partials_.reserve(32);
        if(value != 0.0) {
            partials_.push_back(value);
        }
    }

    //TODO: Adapt this for compilers other than GCC.
    // VS has a pragma
    // Clang seems to do the right thing
    __attribute__((__optimize__("no-associative-math")))
    inline ExactSum& operator()(double x) {
        typedef std::vector<double>::size_type size_type;
        size_type i=0;
        double xsave = x;
        for(size_type j=0;j<partials_.size();++j) {
            double y = partials_[j], hi = x+y;
            double lo = (fabs(x) >= fabs(y)) ?
                y - (double)(hi-x) : x - (double)(hi-y);
            if(lo != 0.0) {
                partials_[i++] = lo;
            }
            x = hi;
        }
        if(x == 0.0) {
            partials_.resize(i);
        } else if(std::isfinite(x)) {
            partials_.resize(i+1);
            partials_[i] = x;
        } else if(std::isfinite(xsave)) {
            // Intermediate overflow has occurred
            failed_ = true;
            partials_.clear();
        } else {
            if(std::isinf(xsave)) {
                inf_sum_ += xsave;
            }
            special_sum_ += xsave;
            partials_.clear();
        }
        return *this;
    }

    __attribute__((__optimize__("no-associative-math")))
    inline double result() const {
        typedef std::vector<double>::size_type size_type;
        if(failed_) {
            return HUGE_VAL;
        } else if(special_sum_ != 0.0) {
            return std::isnan(inf_sum_) ? inf_sum_ : special_sum_; 
        } else if(partials_.empty()) {
            return 0.0;
        }
        size_type j = partials_.size()-1;
        double hi = partials_[j];
        double lo;
        while(j > 0) {
            double x = hi;
            double y = partials_[--j];
            assert(fabs(y) < fabs(x));
            hi = x + y;
            lo = y - (double)(hi - x);
            if(lo != 0.0)
                break;
        }
        if(j > 0 && ((lo < 0.0 && partials_[j-1] < 0.0 ) ||
                     (lo > 0.0 && partials_[j-1] > 0.0))) {
            double y = lo * 2.0;
            double x = hi + y;
            if( y == x-hi ) {
                hi = x;
            }
        }
        return hi;
    }

    inline ExactSum& operator+=(double value) {
        return operator()(value);
    }

    inline ExactSum& add(double value) {
        return operator()(value);
    }

    inline operator double() const {
        return result();
    }
    inline bool operator==(double x) const {
        return (result() == x);
    }

    inline bool failed() const {
        return failed_;
    }
    inline operator bool() const {
        return failed_;
    }
    inline bool operator==(bool b) const {
        return failed_ == b;
    }

private:
    std::vector<double> partials_;
    double special_sum_{0.0};
    double inf_sum_{0.0};
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
