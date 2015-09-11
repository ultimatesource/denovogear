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


#include <dng/matrix.h>
#include <dng/likelihood.h>

#ifndef DNG_LIKLIHOOD_PHI_MIN
#   define DNG_LIKLIHOOD_PHI_MIN 1e-6
#endif

dng::genotype::DirichletMultinomialMixture
::DirichletMultinomialMixture(params_t model_a, params_t model_b) :
    cache_(5), alphas_(5) {
    double a, u, e, m, h;

    f1_ = log(model_a.pi) - log(model_a.pi + model_b.pi);
    f2_ = log(model_b.pi) - log(model_a.pi + model_b.pi);

    // model a
    // Don't allow phi to be zero until we have rewritten the likelihood function
    if(model_a.phi < DNG_LIKLIHOOD_PHI_MIN) {
        model_a.phi = DNG_LIKLIHOOD_PHI_MIN;
    }
    a = (1.0 - model_a.phi) / model_a.phi;
    u = model_a.omega;
    e = (model_a.epsilon * u) / 3.0 / (1.0 - model_a.epsilon * (1.0 - u));

    m = 1.0 - 3.0 * e; // prob of a read that matches a homozygote
    h = (1.0 - 2.0 * e) / 2.0; // prob of a read that matches a heterozygote
    for(int r = 0; r < 5; ++r) {
        for(int g = 0; g < 10; ++g) {
            double tmp[5] = {e, e, e, e, e};
            int g1 = folded_diploid_nucleotides[g][0];
            int g2 = folded_diploid_nucleotides[g][1];
            if(g1 == g2) {
                tmp[g1] = m;
            } else {
                tmp[g1] = h;
                tmp[g2] = h;
            }
            tmp[r] *= u;
            tmp[4] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
            double aa = a / tmp[4];
            for(int x = 0; x < 5; ++x) {
                alphas_[r][g][x][0] = aa * tmp[x];
                alphas_[r][g][x][1] = lgamma(aa * tmp[x]);
            }
        }
    }

    // model b
    if(model_b.phi < DNG_LIKLIHOOD_PHI_MIN) {
        model_b.phi = DNG_LIKLIHOOD_PHI_MIN;
    }
    a = (1.0 - model_b.phi) / model_b.phi;
    u = model_b.omega;
    e = (model_b.epsilon * u) / 3.0 / (1.0 - model_b.epsilon * (1.0 - u));

    m = 1.0 - 3.0 * e; // prob of a read that matches a homozygote
    h = (1.0 - 2.0 * e) / 2.0; // prob of a read that matches a heterozygote
    for(int r = 0; r < 5; ++r) {
        for(int g = 0; g < 10; ++g) {
            double tmp[5] = {e, e, e, e, e};
            int g1 = folded_diploid_nucleotides[g][0];
            int g2 = folded_diploid_nucleotides[g][1];
            if(g1 == g2) {
                tmp[g1] = m;
            } else {
                tmp[g1] = h;
                tmp[g2] = h;
            }
            tmp[r] *= u;
            tmp[4] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
            double aa = a / tmp[4];
            for(int x = 0; x < 5; ++x) {
                alphas_[r][g][x][2] = aa * tmp[x];
                alphas_[r][g][x][3] = lgamma(aa * tmp[x]);
            }
        }
    }

    // construct cache
    for(int r = 0; r < 5; ++r) { // reference
        for(int g = 0; g < 10; ++g) { // genotype
            for(int x = 0; x < 5; ++x) { // nucleotide
                double t1 = 0.0, t2 = 0.0;
                for(int k = 0; k < kCacheSize; ++k) {
                    cache_[r][g][x][2 * k] = t1;
                    cache_[r][g][x][2 * k + 1] = t2;
                    t1 += log(alphas_[r][g][x][0] + k);
                    t2 += log(alphas_[r][g][x][2] + k);
                }
            }
        }
    }
}
