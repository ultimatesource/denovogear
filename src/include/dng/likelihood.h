/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
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
#ifndef DNG_LIKELIHOOD_H
#define DNG_LIKELIHOOD_H

#include <array>
#include <cmath>
#include <memory>
#include <iostream>

#include <boost/math/special_functions/lanczos.hpp>

#include <dng/matrix.h>
#include <dng/utility.h>
#include <dng/depths.h>

namespace dng {
namespace genotype {

namespace detail {
class log_pochhammer {
    typedef boost::math::lanczos::lanczos13m53 Lanczos;
public:
    static constexpr int kCacheSize{512};

    log_pochhammer() { }

    log_pochhammer(double a) : cache_(kCacheSize) {
        // store these values for later usage
        a_ = a;
        ah_ = a - 0.5;
        loga_ = log(a);
        agh_ = a + Lanczos::g() - 0.5;
        la_ = log(Lanczos::lanczos_sum(a));
        logam1_ = loga_-1.0;
        for(int i=0;i<kCacheSize;++i) {
            cache_[i] = lnpoch_pos(i+1);
        }
    }

    double operator()(int n) const {
        assert(n >= 0);
        assert(cache_.size() == kCacheSize); // Check for proper initialization
        if(n <= 0) {
            return 0.0;
        } else if(n <= kCacheSize) {
            return cache_[n-1];
        }
        return lnpoch_pos(n);
    }
private:
    // store these values for later usage
    double a_, loga_, agh_, ah_, la_, logam1_;

    std::vector<double> cache_;

    double lnpoch_pos(double n) const {
        assert(n >= 1.0);
        // if n is small relative to a, we can use a Sterling-derived approximation
        if(n < a_*sqrt(DBL_EPSILON)) {
            return n*logam1_ + (n+ah_)*log1p(n/a_);
        }
        // Use Boost's Lanczos approximation
        double result = ah_*log1p(n/agh_);
        result += log(Lanczos::lanczos_sum(a_+n)) - la_;
        result += n*(log(agh_ + n)-1.0);
        return result;
    }
};

}

class DirichletMultinomialMixture {
public:
    // TODO: Make this configurable with a define
    static const int kCacheSize = 512;

    struct params_t {
        double pi;      // probability of this component
        double phi;     // overdispersion parameter
        double epsilon; // prob of error when homozygote is sequenced
        double omega;   // bias towards reference when heterozygote is sequenced
        // fraction of reads that are ref is omega/(1+omega)
        // Construct a params_t from a string of comma separated values
        params_t(const std::string &str) {
            auto f = utility::parse_double_list(str, ',', 4);
            if(!f.second) {
                throw std::runtime_error("Unable to parse genotype-likelihood parameters. "
                                         "It must be a comma separated list of floating-point numbers.");
            }
            if(f.first.size() != 4) {
                throw std::runtime_error("Wrong number of values for genotype-likelihood parameters. "
                                         "Expected 4; found " + std::to_string(f.first.size()) + ".");
            }
            pi = f.first[0];
            phi = f.first[1];
            epsilon = f.first[2];
            omega = f.first[3];
        }
    };

    std::pair<GenotypeArray, double> operator()(depth_t d, int ref_allele) const;
    
    double operator()(const pileup::AlleleDepths& depths, const std::vector<size_t> &indexes, IndividualVector::iterator output) const;

    DirichletMultinomialMixture(params_t model_a, params_t model_b);

protected:

    // NOTE: a = reference; b = genotype; c = nucleotide; d = depth
    // NOTE: cache_[a][b][c].first(d)  = sum alpha1[a][b][c]+x for x in [0,d)
    // NOTE: cache_[a][b][c].second(d) = sum alpha2[a][b][c]+x for x in [0,d)
    typedef detail::log_pochhammer cache_type;
    typedef std::vector<std::array<std::array<std::pair<cache_type,cache_type>, 5>, 10>> cache_t;

    cache_t cache_;
    double f1_, f2_; // log(f) and log(1-f)
};

} // namespace genotype
} // namespace dng

#endif // DNG_LIKELIHOOD_H

