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


#include <dng/matrix.h>
#include <dng/likelihood.h>

//  The smallest positive value such that (1-phi)/phi + 1 != (1-phi)/phi
#ifndef DNG_LIKLIHOOD_PHI_MIN
#   define DNG_LIKLIHOOD_PHI_MIN (DBL_EPSILON/2.0)
#endif

namespace dng {
namespace genotype {

DirichletMultinomialMixture::DirichletMultinomialMixture(params_t model_a, params_t model_b) :
    cache_(5) {
    double a, u, e, m, h;

    assert(0.0 < model_a.pi && 0.0 < model_b.pi);
    assert(0.0 <= model_a.phi && model_a.phi <= 1.0 );
    assert(0.0 <= model_b.phi && model_b.phi <= 1.0 );
    assert(0.0 <= model_a.epsilon && model_a.epsilon <= 1.0 );
    assert(0.0 <= model_b.epsilon && model_b.epsilon <= 1.0 );
    assert(0.0 < model_a.omega && 0.0 < model_b.omega);

    // Calculate log(mixing proportions) and ensure that they are normalized
    f1_ = log(model_a.pi) - log(model_a.pi + model_b.pi);
    f2_ = log(model_b.pi) - log(model_a.pi + model_b.pi);

    // model a
    // Don't allow phi to be zero until we have rewritten the likelihood function
    if(model_a.phi < DNG_LIKLIHOOD_PHI_MIN) {
        model_a.phi = DNG_LIKLIHOOD_PHI_MIN;
    }
    a = (1.0 - model_a.phi) / model_a.phi;
    u = model_a.omega;
    e = model_a.epsilon/3.0;

    m = 1.0 - 3.0 * e; // prob of a read that matches a homozygote
    h = 1.0 - 2.0 * e; // prob of a read that matches a heterozygote
    for(int r = 0; r < 5; ++r) {
        for(int g = 0; g < 10; ++g) {
            double tmp[5] = {e, e, e, e, e};
            int g1 = folded_diploid_nucleotides[g][0];
            int g2 = folded_diploid_nucleotides[g][1];
            if(g1 == g2) {
                tmp[g1] = m;
            } else if(g1 == r) {
                tmp[g1] = h*u/(1.0+u);
                tmp[g2] = h*1.0/(1.0+u);
            } else if(g2 == r) {
                tmp[g1] = h*1.0/(1.0+u);
                tmp[g2] = h*u/(1.0+u);
            } else {
                tmp[g1] = h/2.0;
                tmp[g2] = h/2.0;
            }
            tmp[4] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
            double aa = a / tmp[4];
            for(int x = 0; x < 5; ++x) {
                // TODO: Is it worth eliminating redundant calculations?
                cache_[r][g][x].first = cache_type{aa*tmp[x]};
            }
        }
    }

    // model b
    if(model_b.phi < DNG_LIKLIHOOD_PHI_MIN) {
        model_b.phi = DNG_LIKLIHOOD_PHI_MIN;
    }
    a = (1.0 - model_b.phi) / model_b.phi;
    u = model_b.omega;
    e = model_b.epsilon/3.0;

    m = 1.0 - 3.0 * e; // prob of a read that matches a homozygote
    h = 1.0 - 2.0 * e; // prob of a read that matches a heterozygote
    for(int r = 0; r < 5; ++r) {
        for(int g = 0; g < 10; ++g) {
            double tmp[5] = {e, e, e, e, e};
            int g1 = folded_diploid_nucleotides[g][0];
            int g2 = folded_diploid_nucleotides[g][1];
            if(g1 == g2) {
                tmp[g1] = m;
            } else if(g1 == r) {
                tmp[g1] = h*u/(1.0+u);
                tmp[g2] = h*1.0/(1.0+u);
            } else if(g2 == r) {
                tmp[g1] = h*1.0/(1.0+u);
                tmp[g2] = h*u/(1.0+u);
            } else {
                tmp[g1] = h/2.0;
                tmp[g2] = h/2.0;
            }
            tmp[4] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
            double aa = a / tmp[4];
            for(int x = 0; x < 5; ++x) {
                // TODO: Is it worth eliminating redundant calculations?
                cache_[r][g][x].second = cache_type{aa*tmp[x]};
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

double DirichletMultinomialMixture::operator()(
    const pileup::AlleleDepths& depths, const std::vector<size_t> &indexes,
    IndividualVector::iterator output) const
{
    const auto last = output + indexes.size();

    int ref_index = depths.type_info().reference;
    int width = depths.type_info().width;
    int gt_width = depths.type_gt_info().width;
    // resize output to hold genotypes
    for(auto first = output; first != last; ++first ) {
        first->resize(gt_width);
    }

    // calculate log_likelihoods for each genotype
    for(int g=0;g<gt_width;++g) {
        int real_g = depths.type_gt_info().indexes[g];
        auto &cache = cache_[ref_index][real_g];
        for(size_t u = 0; u < indexes.size(); ++u) {
            size_t lib = indexes[u];
            // create a reference to the output site
            auto &out = *(output+u);
            double lh1 = f1_, lh2 = f2_;
            int read_count = 0;
            for(int i=0;i<width;++i) {
                // get the depth of nucleotide i from library lib
                int d = depths(i,lib);
                read_count += d;
                // find which nucleotide this depth refers to
                int j = depths.type_info().indexes[i];
                lh1 += cache[j].first(d);
                lh2 += cache[j].second(d);
            }
            lh1 -= cache[4].first(read_count);
            lh2 -= cache[4].second(read_count);
            out[g] = log_sum(lh1,lh2);
            assert(std::isfinite(out[g]));
        }
    }
    // rescale output to hold genotypes
    double scale = 0.0;
    for(auto first = output; first != last; ++first ) {
        double scaleg = first->maxCoeff();
        scale += scaleg;
        *first = (*first - scaleg).exp();
    }
    return scale;
}

std::pair<GenotypeArray, double> DirichletMultinomialMixture::operator()(depth_t d, int ref_allele) const {
    GenotypeArray log_ret{10};
    int read_count = d.counts[0] + d.counts[1] +
                     d.counts[2] + d.counts[3];
    for(int i = 0; i < 10; ++i) {
        auto &cache = cache_[ref_allele][i];
        double lh1 = f1_, lh2 = f2_;
        for(int j=0;j<4;++j) {
            lh1 += cache[j].first(d.counts[j]);
            lh2 += cache[j].second(d.counts[j]);
        }
        lh1 -= cache[4].first(read_count);
        lh2 -= cache[4].second(read_count);

        log_ret[i] = log_sum(lh1,lh2);
    }
    double scale = log_ret.maxCoeff();
    return {(log_ret - scale).exp(), scale};
    }

}} // namespace dng::genotype
