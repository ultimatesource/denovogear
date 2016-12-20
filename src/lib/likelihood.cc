/*
 * Copyright (c) 2014-2016 Reed A. Cartwright
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


#include <dng/likelihood.h>

//  The smallest positive value such that (1-phi)/phi + 1 != (1-phi)/phi
#ifndef DNG_LIKLIHOOD_PHI_MIN
#   define DNG_LIKLIHOOD_PHI_MIN (DBL_EPSILON/2.0)
#endif

namespace dng {
namespace genotype {

std::array<double,5> make_alphas(int reference, int genotype, double phi, double epsilon, double omega) {
    assert(0 <= reference && reference <= 4);
    assert(0 <= genotype && genotype <= 9);
    assert(0.0 <= phi && phi <= 1.0 );
    assert(0.0 <= epsilon && epsilon <= 1.0 );
    assert(0.0 <= omega);

    // Don't allow phi to be zero until we have rewritten the likelihood function
    if(phi < DNG_LIKLIHOOD_PHI_MIN) {
        phi = DNG_LIKLIHOOD_PHI_MIN;
    }
    double a = (1.0 - phi) / phi;
    double u = omega;
    double e = epsilon/3.0;
    double m = 1.0 - 3.0 * e; // prob of a read that matches a homozygote
    double h = 1.0 - 2.0 * e; // prob of a read that matches a heterozygote

    std::array<double,5> ret = {e,e,e,e,e};
    int g1 = folded_diploid_nucleotides[genotype][0];
    int g2 = folded_diploid_nucleotides[genotype][1];
    if(g1 == g2) {
        ret[g1] = m;
    } else if(g1 == reference) {
        ret[g1] = h*u/(1.0+u);
        ret[g2] = h*1.0/(1.0+u);
    } else if(g2 == reference) {
        ret[g1] = h*1.0/(1.0+u);
        ret[g2] = h*u/(1.0+u);
    } else {
        ret[g1] = h/2.0;
        ret[g2] = h/2.0;
    }
    double total = ret[0] + ret[1] + ret[2] + ret[3];
    for(int i=0;i<4;++i) {
        ret[i] = ret[i]/total * a;
    }
    ret[4] = a;
    return ret;
}

DirichletMultinomialMixture::DirichletMultinomialMixture(params_t model_a, params_t model_b) {
    assert(0.0 < model_a.pi && 0.0 < model_b.pi);
    // Calculate log(mixing proportions) and ensure that they are normalized
    f1_ = log(model_a.pi) - log(model_a.pi + model_b.pi);
    f2_ = log(model_b.pi) - log(model_a.pi + model_b.pi);

    // Construct our log_pochhammer functors and cache some of their outputs
    for(int reference=0; reference<5; ++reference) {
        for(int genotype=0; genotype<10; ++genotype) {
            auto alphas0 = make_alphas(reference, genotype,
                model_a.phi, model_a.epsilon, model_a.omega);
            auto alphas1 = make_alphas(reference, genotype,
                model_b.phi, model_b.epsilon, model_b.omega);            
            for(int nucleotide=0; nucleotide<5; ++nucleotide) {
                // Construct functors for log_gamma(a+n)-log_gamma(a)
                // where a = alpha(reference, genotype, nucleotide)
                // where n = depth of that nucleotide
                // nucleotide = 4 refers to the total coverage
                // first component
                models_[reference][nucleotide][0][genotype] = detail::log_pochhammer{alphas0[nucleotide]};
                // second component
                models_[reference][nucleotide][1][genotype] = detail::log_pochhammer{alphas1[nucleotide]};
                
                // Cache the results of the functor for n = [0,kCacheSize)
                for(int i=0;i<cache_.size();++i) {
                    if(nucleotide < 4) {
                        cache_[i][reference][nucleotide][0][genotype] =
                            models_[reference][nucleotide][0][genotype](i);
                        cache_[i][reference][nucleotide][1][genotype] =
                            models_[reference][nucleotide][1][genotype](i);
                    } else {
                        cache_[i][reference][nucleotide][0][genotype] =
                            -models_[reference][nucleotide][0][genotype](i);
                        cache_[i][reference][nucleotide][1][genotype] =
                            -models_[reference][nucleotide][1][genotype](i);                        
                    }
                }
            }
        }
    }
}

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

std::pair<GenotypeArray, double> DirichletMultinomialMixture::operator()(
        depth_t d, int ref_allele, int ploidy) const {

    // how many genotypes will we calculate based on ploidy?
    const int sz = (ploidy==2) ? 10 : 4;
    
    // Our return variable
    GenotypeArray log_ret{sz};

    const int read_count = d.counts[0] + d.counts[1] +
                     d.counts[2] + d.counts[3];

    // use a single temporary buffer for speed
    cache_type temp[20];
    std::fill(&temp[0],&temp[sz],f1_);
    std::fill(&temp[10],&temp[10+sz],f2_);

    // Loop through the depths and add log_gamma(a+n)-log_gamma(a) to each temporary genotype
    // The order that this loop proceeds has been chosen to hopefully
    // utilize cache/memory effectively.
    for(int nucleotide = 0; nucleotide < 5; ++nucleotide) {
        int count = (nucleotide < 4) ? d.counts[nucleotide] : read_count;
        if(count < cache_.size()) {
            // If count is reasonable, use cached constants
            const auto &cache = cache_[count][ref_allele][nucleotide];
            for(int genotype = 0; genotype < sz; ++genotype) {
                temp[genotype] += cache[0][genotype];
            }
            for(int genotype = 0; genotype < sz; ++genotype) {
                temp[10+genotype] += cache[1][genotype];
            }
        } else if(nucleotide < 4) {
            // If count is big, use our model functor
            const auto &model = models_[ref_allele][nucleotide];
            for(int genotype = 0; genotype < sz; ++genotype) {
                temp[genotype] += model[0][genotype](count);
            }
            for(int genotype = 0; genotype < sz; ++genotype) {
                temp[10+genotype] += model[1][genotype](count);
            }
        } else {
            // If nucleotide == 4, we have to use a subtraction
            const auto &model = models_[ref_allele][nucleotide];
            for(int genotype = 0; genotype < sz; ++genotype) {
                temp[genotype] -= model[0][genotype](count);
            }
            for(int genotype = 0; genotype < sz; ++genotype) {
                temp[10+genotype] -= model[1][genotype](count);
            }
        }
    }
    // Calculate log_likelihoods for genotypes
    for(int genotype = 0; genotype < sz; ++genotype) {
        log_ret[genotype] = log_sum_exact(temp[genotype],temp[10+genotype]);
        assert(std::isfinite(log_ret[genotype]));
    }
    // Scale and calculate likelihoods
    double scale = log_ret.maxCoeff();
    return {(log_ret - scale).exp(), scale};
}

std::pair<GenotypeArray, double> DirichletMultinomialMixture::operator()(
        const pileup::AlleleDepths& depths, size_t pos, int ploidy) const {
    int ref_allele = depths.type_info().reference;
    int width = depths.type_info().width;
    int sz = (ploidy == 2) ? depths.type_gt_info().width : depths.type_info().width;
    const char *indexes = (ploidy == 2) ? &depths.type_gt_info().indexes[0] : &depths.type_info().indexes[0];

    // Our return variable
    GenotypeArray log_ret{sz};

    // use a single temporary buffer for speed
    cache_type temp[20];
    std::fill(&temp[0],&temp[sz],f1_);
    std::fill(&temp[10],&temp[10+sz],f2_);

    int total = 0;
    for(int i = 0; i < width; ++i) {
        const int nucleotide = depths.type_info().indexes[i];
        const int count = depths(pos,i);
        total += count;
        if(count < cache_.size()) {
            // If count is reasonable, use cached constants
            const auto &cache = cache_[count][ref_allele][nucleotide];
            for(int j = 0; j < sz; ++j) {
                int genotype = indexes[j];
                temp[j] += cache[0][genotype];
            }
            for(int j = 0; j < sz; ++j) {
                int genotype = indexes[j];
                temp[10+j] += cache[1][genotype];
            }
        } else {
            const auto &model = models_[ref_allele][nucleotide];
            for(int j = 0; j < sz; ++j) {
                int genotype = indexes[j];
                temp[j] += model[0][genotype](count);
            }
            for(int j = 0; j < sz; ++j) {
                int genotype = indexes[j];
                temp[10+j] += model[1][genotype](count);
            }            
        }
    }
    const int count = total;
    const int nucleotide = 4;
    if(count < cache_.size()) {
        // If count is reasonable, use cached constants
        const auto &cache = cache_[count][ref_allele][nucleotide];
        for(int j = 0; j < sz; ++j) {
            int genotype = indexes[j];
            temp[j] += cache[0][genotype];
        }
        for(int j = 0; j < sz; ++j) {
            int genotype = indexes[j];
            temp[10+j] += cache[1][genotype];
        }
    } else {
        const auto &model = models_[ref_allele][nucleotide];
        for(int j = 0; j < sz; ++j) {
            int genotype = indexes[j];
            temp[j] -= model[0][genotype](count);
        }
        for(int j = 0; j < sz; ++j) {
            int genotype = indexes[j];
            temp[10+j] -= model[1][genotype](count);
        }            
    }
    // Calculate log_likelihoods for genotypes
    for(int j = 0; j < sz; ++j) {
        log_ret[j] = log_sum_exact(temp[j],temp[10+j]);
        assert(std::isfinite(log_ret[j]));
    }
    // Scale and calculate likelihoods
    double scale = log_ret.maxCoeff();
    return {(log_ret - scale).exp(), scale};
}


} // namespace genotype
} // namespace dng
