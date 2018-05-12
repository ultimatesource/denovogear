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

enum struct Mode {
    Likelihood = 0,
    LogLikelihood
};

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

std::array<double,8> make_alphas(double over_dispersion_hom, 
        double over_dispersion_het, double ref_bias,
        double error_rate, double error_alleles);

} // namespace detail

class DirichletMultinomial {
public:
    using depths_const_reference_type = dng::pileup::allele_depths_ref_t::const_reference;

    DirichletMultinomial(double over_dispersion_hom, 
        double over_dispersion_het, double sequencing_bias,
        double error_rate, double k_alleles);

    template<typename Range>
    double operator()(const Range& ad, int num_obs_alleles, Mode mode, int ploidy,
        GenotypeArray *output) const;

    double over_dispersion_hom() const { return over_dispersion_hom_; }
    double over_dispersion_het() const { return over_dispersion_het_; }
    double error_rate() const { return error_rate_; }
    double sequencing_bias() const { return sequencing_bias_; }
    double k_alleles() const { return k_alleles_; }

protected:

    template<typename Range>
    void LogHaploidGenotypes(const Range& ad, int num_obs_alleles, GenotypeArray *output) const;

    template<typename Range>
    void LogDiploidGenotypes(const Range& ad, int num_obs_alleles, GenotypeArray *output) const;

    double over_dispersion_hom_; // overdispersion of homozygotes
    double over_dispersion_het_; // overdispersion of heterozygotes
    double error_rate_;      // prob of error when homozygote is sequenced
    double sequencing_bias_;        // bias towards reference when heterozygote is sequenced
    double k_alleles_;

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

    friend std::array<double,8> detail::make_alphas(double over_dispersion_hom, 
        double over_dispersion_het, double ref_bias,
        double error_rate, double k_alleles);
};

template<typename Range>
void DirichletMultinomial::LogDiploidGenotypes(const Range& ad, int num_obs_alleles, GenotypeArray *output) const {    
    assert(num_obs_alleles >= 1);
    assert(output != nullptr);
    const int sz = num_obs_alleles*(num_obs_alleles+1)/2;

    auto &ret = *output;
    ret.setZero(sz);

    int total = 0, pos=0;
    for(auto d : ad) {
        assert( d >= 0 );
        total += d;

        for(int a=0,gt=0; a < num_obs_alleles; ++a) {
            for(int b=0; b < a; ++b) {
                // Heterozygotes
                if(b == 0) {
                    // gt = 0/a
                    if(pos == b) {
                        ret[gt++] += pochhammer(alpha::HET_REF, d);    
                    } else if(pos == a) {
                        ret[gt++] += pochhammer(alpha::HET_ALT, d);
                    } else {
                        ret[gt++] += pochhammer(alpha::HET_ERROR, d);
                    }
                } else {
                    // gt = b/a
                    if(pos == b || pos == a) {
                        ret[gt++] += pochhammer(alpha::HET_ALTALT, d); 
                    } else {
                        ret[gt++] += pochhammer(alpha::HET_ERROR, d);
                    }
                }
            }
            // Homozygotes
            if(pos == a) {
                ret[gt++] += pochhammer(alpha::HOM_MATCH, d);
            } else {
                ret[gt++] += pochhammer(alpha::HOM_ERROR, d);
            }
        }
        ++pos;
    } 

    double hom_total = pochhammer(alpha::HOM_TOTAL, total);
    ret[0] -= hom_total;
    if(num_obs_alleles > 1) {
        double het_total = pochhammer(alpha::HET_TOTAL, total);
        for(int a=1,gt=1; a < num_obs_alleles; ++a) {
            for(int b=0; b < a; ++b) {
                ret[gt++] -= het_total;
            }
            ret[gt++] -= hom_total;
        }
    }
}

template<typename Range>
void DirichletMultinomial::LogHaploidGenotypes(const Range& ad, int num_obs_alleles, GenotypeArray *output) const {
    assert(num_obs_alleles >= 1); 
    assert(output != nullptr);
    const int sz = num_obs_alleles;

    auto &ret = *output;
    ret.setZero(sz);

    int total = 0, pos = 0;
    for(auto d : ad) {
        assert( d >= 0 );
        total += d;
        for(int gt=0; gt < sz; ++gt) {
            ret[gt] += pochhammer((gt == pos) ? alpha::HOM_MATCH : alpha::HOM_ERROR, d) ;
        }
        ++pos;
    }
    double hom_total = pochhammer(alpha::HOM_TOTAL, total);
    for(int gt=0; gt < sz; ++gt) {
        ret[gt] -=  hom_total;
    }
}


template<typename Range>
double DirichletMultinomial::operator()(const Range &ad, int num_obs_alleles, Mode mode, int ploidy,
    GenotypeArray *output) const
{
    assert(output != nullptr); 
    assert(num_obs_alleles >= 1);
    assert(ploidy == 1 || ploidy == 2);

    if(ploidy == 2) {
        LogDiploidGenotypes(ad,num_obs_alleles,output);
    } else {
        LogHaploidGenotypes(ad,num_obs_alleles,output);
    }

    // Scale and calculate likelihoods
    double scale = output->maxCoeff();
    if(mode == Mode::Likelihood) {
        *output = (*output - scale).exp();        
    } else {
        *output = (*output - scale);
    }
    return scale;
}


} // namespace genotype

using Genotyper = genotype::DirichletMultinomial;

} // namespace dng

#endif // DNG_LIKELIHOOD_H
