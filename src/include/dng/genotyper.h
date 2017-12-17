/*
 * Copyright (c) 2014-2017 Reed A. Cartwright
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
#ifndef DNG_GENOTYPER_H
#define DNG_GENOTYPER_H

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
    log_pochhammer() { }

    log_pochhammer(double a) {
        // store these values for later usage
        assert(a > 0.0);
        a_ = a;
        ah_ = a - 0.5;
        agh_ = a + Lanczos::g() - 0.5;
        la_ = log(Lanczos::lanczos_sum(a));
        logam1_ = log(a)-1.0;
    }

    double operator()(int n) const {
        return lnpoch_pos(n);
    }
private:
    // store these values for later usage
    double a_, agh_, ah_, la_, logam1_;

    double lnpoch_pos(double n) const {
        assert(n >= 0.0);
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
} // namespace detail

class DirichletMultinomial {
public:
    using allele_depths_t = dng::pileup::allele_depths_t;

    DirichletMultinomial(double over_dispersion_hom, 
        double over_dispersion_het, double ref_bias,
        double error_rate, double error_entropy);

    std::pair<GenotypeArray, double> operator()(
        allele_depths_t::const_reference ad, int num_alts, int ploidy=2) const;

    double over_dispersion_hom() const { return over_dispersion_hom_; }
    double over_dispersion_het() const { return over_dispersion_het_; }
    double error_rate() const { return error_rate_; }
    double ref_bias() const { return ref_bias_; }

protected:
    GenotypeArray LogHaploidGenotypes(allele_depths_t::const_reference ad, int num_alts) const;
    GenotypeArray LogDiploidGenotypes(allele_depths_t::const_reference ad, int num_alts) const;

    double over_dispersion_hom_; // overdispersion of homozygotes
    double over_dispersion_het_; // overdispersion of heterozygotes
    double error_rate_;      // prob of error when homozygote is sequenced
    double ref_bias_;        // bias towards reference when heterozygote is sequenced
    double error_entropy_;

    using cache_t = std::vector<std::array<double,8>>;
    using pochhammers_t = std::array<detail::log_pochhammer,8>; 
    static constexpr int CACHE_SIZE = 512;
    cache_t cache_{CACHE_SIZE};
    pochhammers_t pochhammers_;

    enum struct alpha {
        HOM_MATCH = 0, HOM_ERROR,
        HET_REF, HET_ALT, HET_ERROR,
        HET_ALTALT, HOM_TOTAL, HET_TOTAL
    };

    inline double pochhammer(alpha a, int n) const {
        assert(n >= 0);
        int t = static_cast<int>(a);
        assert(0 <= t && t < 8);
        return (n < CACHE_SIZE) ? cache_[n][t] : pochhammers_[t](n);
    }

    friend std::array<double,8> make_alphas(double over_dispersion_hom, 
        double over_dispersion_het, double ref_bias,
        double error_rate, double error_entropy);
};

} // namespace genotype

using Genotyper = genotype::DirichletMultinomial;

} // namespace dng

#endif // DNG_LIKELIHOOD_H
