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

/* K-Alleles Genotyping Model

# Model Parameters
 - e: error rate
 - w: ref bias in 0/1 heterozygotes
 - o1: overdispersion in homozygotes  
 - o2: overdispersion in heterozygotes 

# Intermediate Parameters
 - p1 = (1-e) --- the proportion of base-calls that match a homozygote
 - p2 = e/(k-1) --- the probability of a specific erroneous base call
 - p3 = (1-e+e/(k-1)) -- the proportion of base-calls that match a heterozyote
 - q = w/(1+w) --- the fraction of matching reads in a ref/alt het that match ref

# Homozygote Alphas
  - o1*p1 # matches
  - o1*p2 # mismatches

# Heterozygote Alphas
  - o2*p3*q # ref matches
  - o2*p3*(1-q) # alt matches
  - o2*p3*0.5 # matches for a/b hets
  - o2*p2 # mismatches
*/


#include <dng/genotyper.h>

// The smallest positive value such that (1-phi)/phi + 1 != (1-phi)/phi
#ifndef DNG_LIKLIHOOD_PHI_MIN
#   define DNG_LIKLIHOOD_PHI_MIN (DBL_EPSILON/2.0)
#endif
// Use ups for lowest error rate, which prevents alpha=0 for log_pochhammer
#ifndef DNG_LIKLIHOOD_EPSILON_MIN
#   define DNG_LIKLIHOOD_EPSILON_MIN (DBL_EPSILON/2.0)
#endif


namespace dng {
namespace genotype {

namespace detail {
// lookup table is stored in another file
#include "log1p_exp_neg_table.tcc"

// calculate m+log(1+exp(-x))
inline double log1p_exp_neg(double x)  {
    assert(x >= 0.0);
    if(x < 2.0) {
        // use a Taylor series of degree 8 at x=0
        double xx = x*x;
        return M_LN2-x/2.0 + (xx*(0.125+ xx*(-1.0/192.0+ (1.0/2880.0 - 17.0/645120.0*xx)*xx)));
    } else if(x < 7.0) {
        // use a Taylor series of degree 8 at x=4.694144053627953
        return 0.7660035463487451 + x*(-0.6446774487882012 + x*(0.24816087438710718 + 
            x*(-0.05625031169464227 + x*(0.008050735924081637 + x*(-0.0007225510231045971 +
            (0.00003741637881581404 - 8.576355342428642e-7*x)*x)))));
    } else if(x < log1p_exp_neg_table_end) {
        x = x-7.0;
        int n = static_cast<int>(x/log1p_exp_neg_table_step);
        assert(n < 1023);
        double f0 = log1p_exp_neg_table[n];
        double f1 = log1p_exp_neg_table[n+1];
        double dx = x-(n*log1p_exp_neg_table_step);
        double df = f1-f0;
        return f0 + df*dx/log1p_exp_neg_table_step;
    }
    return 0.0;
}

inline double log_sum(double a, double b) {
     double x = fabs(a-b);
     double m = std::max(a,b);
     return m + log1p_exp_neg(x);
}

inline double log_sum_exact(double a, double b) {
    return log1p(exp(-fabs(a-b))) + std::max(a,b);
}
}

std::array<double,8> detail::make_alphas(double over_dispersion_hom, 
        double over_dispersion_het, double ref_bias,
        double error_rate, double error_entropy) {
    assert(0.0 <= over_dispersion_hom && over_dispersion_hom <= 1.0 );
    assert(0.0 <= over_dispersion_het && over_dispersion_het <= 1.0 );
    assert(0.0 <= ref_bias );
    assert(0.0 <= error_rate && error_rate <= 1.0 );
    assert(0.0 <= error_entropy );

    // Adjust parameters to prevent zeros, which can cause issues with the 
    // current implementation
    over_dispersion_hom = std::max(over_dispersion_hom, DNG_LIKLIHOOD_PHI_MIN);
    over_dispersion_het = std::max(over_dispersion_het, DNG_LIKLIHOOD_PHI_MIN);
    error_rate = std::max(error_rate, DNG_LIKLIHOOD_EPSILON_MIN);

    std::array<double,8> ret;

    // assume three possible errors
    double p1 = 1.0-error_rate;
    double p2 = error_rate/exp(error_entropy);
    double p3 = p1+p2;
    double q = ref_bias/(1.0+ref_bias);

    using alpha = DirichletMultinomial::alpha;

    double a1 = (1.0-over_dispersion_hom)/over_dispersion_hom;
    double a2 = (1.0-over_dispersion_het)/over_dispersion_het;

    ret[(int)alpha::HOM_MATCH] = a1*p1;
    ret[(int)alpha::HOM_ERROR] = a1*p2;
    
    ret[(int)alpha::HET_REF] = a2*p3*q;
    ret[(int)alpha::HET_ALT] = a2*p3*(1.0-q);
    ret[(int)alpha::HET_ALTALT] = a2*p3*0.5;
    ret[(int)alpha::HET_ERROR] = a2*p2;

    ret[(int)alpha::HOM_TOTAL] = a1;
    ret[(int)alpha::HET_TOTAL] = a2;

    return ret;
}

using detail::log_sum;

DirichletMultinomial::DirichletMultinomial(double over_dispersion_hom, 
        double over_dispersion_het, double ref_bias, double error_rate, double error_entropy) : 
    over_dispersion_hom_{over_dispersion_hom},
    over_dispersion_het_{over_dispersion_het},
    ref_bias_{ref_bias}, error_rate_{error_rate},
    error_entropy_{error_entropy}
{
    auto alphas = detail::make_alphas(over_dispersion_hom, over_dispersion_hom, 
        ref_bias, error_rate, error_entropy);

    for(int i=0; i<pochhammers_.size(); ++i) {
        pochhammers_[i] = pochhammers_t::value_type{alphas[i]};
        for(int n=0; n<CACHE_SIZE; ++n) {
            cache_[n][i] = pochhammers_[i](n);
        }
    }
}

GenotypeArray DirichletMultinomial::LogDiploidGenotypes(depths_const_reference_type ad, int num_alts) const
{    
    assert(num_alts >= 0);
    const int num_alleles = num_alts + 1;
    const int sz = num_alleles*(num_alleles+1)/2;

    GenotypeArray ret{sz};
    ret.setZero();

    int total = 0;
    for(int pos=0; pos < ad.size(); ++pos) {
        assert( ad[pos] >= 0 );
        int d = ad[pos];
        total += d;

        for(int a=0,gt=0; a < num_alleles; ++a) {
            for(int b=0; b < a; ++b,++gt) {
                // Heterozygotes
                if(b == 0) {
                    // gt = 0/a
                    ret[gt] += pochhammer(
                        (pos == 0) ? alpha::HET_REF :
                        (pos == a) ? alpha::HET_ALT :
                                     alpha::HET_ERROR ,
                        d );
                } else {
                    // gt = b/a
                    ret[gt] += pochhammer(
                        (pos == a || pos == b) ? alpha::HET_ALTALT :
                                                 alpha::HET_ERROR ,
                        d );
                }
            }
            // Homozygotes
            ret[gt++] += pochhammer(
                (a == pos) ? alpha::HOM_MATCH :
                             alpha::HOM_ERROR ,
                d);
        }
    } 

    double het_total = pochhammer(alpha::HET_TOTAL, total);
    double hom_total = pochhammer(alpha::HOM_TOTAL, total);
    for(int a=0,gt=0; a < num_alleles; ++a) {
        for(int b=0; b < a; ++b,++gt) {
            ret[gt] -= het_total;
        }
        ret[gt++] -= hom_total;
    }

    return ret;
}

GenotypeArray DirichletMultinomial::LogHaploidGenotypes(depths_const_reference_type ad, int num_alts) const {
    assert(num_alts >= 0);
    const int num_alleles = num_alts + 1;    
    const int sz = num_alleles;

    GenotypeArray ret{sz};
    ret.setZero();

    int total = 0;
    for(int pos=0; pos < ad.size(); ++pos) {
        assert( ad[pos] >= 0 );
        int d = ad[pos];
        total += d;
        for(int gt=0; gt < sz; ++gt) {
            ret[gt] += pochhammer((gt == pos) ? alpha::HOM_MATCH : alpha::HOM_ERROR, d) ;
        }
    }
    for(int gt=0; gt < sz; ++gt) {
        ret[gt] -= pochhammer(alpha::HOM_TOTAL, total) ;
    }
    return ret;
}


std::pair<GenotypeArray, double> DirichletMultinomial::operator()(
        depths_const_reference_type ad, int num_alts, int ploidy) const
{
    assert(ploidy == 1 || ploidy == 2);

    auto log_ret = (ploidy == 2) ? LogDiploidGenotypes(ad,num_alts)
                                 : LogHaploidGenotypes(ad,num_alts);

    // Scale and calculate likelihoods
    double scale = log_ret.maxCoeff();

    return {(log_ret - scale).exp(), scale};
}



} // namespace genotype
} // namespace dng
