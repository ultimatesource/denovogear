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

#include <dng/matrix.h>
#include <dng/utilities.h>

namespace dng {
namespace genotype {

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
            auto f = util::parse_double_list(str, ',', 4);
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

    std::pair<GenotypeArray, double> operator()(depth_t d, int ref_allele) {
        GenotypeArray log_ret{10};
        int read_count = d.counts[0] + d.counts[1] +
                         d.counts[2] + d.counts[3];
        for(int i = 0; i < 10; ++i) {
            auto &cache = cache_[ref_allele][i];
            double lh1 = f1_, lh2 = f2_;
            for(int j : {0, 1, 2, 3}) {
                if(d.counts[j] < kCacheSize) {
                    // access coefs that are in the cache
                    lh1 += cache[j][2 * d.counts[j]]; // component 1
                    lh2 += cache[j][2 * d.counts[j] + 1]; // component 2
                } else {
                    // if the value is not in the cache, calculate it
                    lh1 += lgamma(alphas_[ref_allele][i][j][0] + d.counts[j])
                           - alphas_[ref_allele][i][j][1];
                    lh2 += lgamma(alphas_[ref_allele][i][j][2] + d.counts[j])
                           - alphas_[ref_allele][i][j][3];
                }
            }
            if(read_count < kCacheSize) {
                lh1 -= cache[4][2 * read_count];
                lh2 -= cache[4][2 * read_count + 1];
            } else {
                lh1 -= lgamma(alphas_[ref_allele][i][4][0] + read_count)
                       - alphas_[ref_allele][i][4][1];
                lh2 -= lgamma(alphas_[ref_allele][i][4][2] + read_count)
                       - alphas_[ref_allele][i][4][3];
            }
            log_ret[i] = (lh2 < lh1) ? lh1 + log1p(exp(lh2 - lh1)) :
                         lh2 + log1p(exp(lh1 - lh2)) ;
        }
        double scale = log_ret.maxCoeff();
        return std::make_pair((log_ret - scale).exp(), scale);
    }

    DirichletMultinomialMixture(params_t model_a, params_t model_b);

protected:
    // NOTE: cache_[a][b][c][d] = sum alpha[a][b][c]+x for x in [0,d/2)
    // NOTE: alpha_[a][b][c] = {alpha1, lgamma(alpha1), alpha2, lgamma(alpha2)}
    typedef std::array<double, 2 * kCacheSize> cache_type;
    typedef std::vector<std::array<std::array<cache_type, 5>, 10>> cache_t;
    typedef std::vector<std::array<std::array<std::array<double, 4>, 5>, 10>>
            alphas_t;

    cache_t cache_;
    alphas_t alphas_;
    double f1_, f2_; // log(f) and log(1-f)
};

} // namespace genotype
} // namespace dng

#endif // DNG_LIKELIHOOD_H

