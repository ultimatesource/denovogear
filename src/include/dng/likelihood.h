/*
 * Copyright (c) 2014 Reed A. Cartwright
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

#include <dng/matrix.h>

namespace dng {
namespace genotype {

class DirichletMultinomialMixture {
public:
	// TODO: Make this configurable with a define
	static const int kCacheSize = 512; 
	
	struct params_t {
		double pi; // probability of this component
		double phi; // overdispersion parameter
		double epsilon; // prob of error when homozygote is sequenced
		double omega; // bias towards reference when heterozygote is sequenced
	};
	
	Vector10d operator()(depth_t d, int ref_allele) {
		Vector10d log_ret{10};
		int read_count = d.counts[0] + d.counts[1] +
		                      d.counts[2] + d.counts[3];
		for(int i=0;i<10;++i) {
			auto &cache = cache_[ref_allele][i];
			double lh1 = f1_, lh2 = f2_;
			for(int j : {0,1,2,3} ) {
				if(d.counts[j] < kCacheSize) {
					// access coefs that is in the cache
					lh1 += cache[j][2*d.counts[j]];   // component 1
					lh2 += cache[j][2*d.counts[j]+1]; // component 2
				} else {
					// if the value is not in the cache, calculate it
					lh1 += lgamma(alphas_[ref_allele][i][j][0]+d.counts[j])
							- alphas_[ref_allele][i][j][1];
					lh2 += lgamma(alphas_[ref_allele][i][j][2]+d.counts[j])
							- alphas_[ref_allele][i][j][3];
				}
			}
			if(read_count < kCacheSize) {
				lh1 -= cache[4][2*read_count];
				lh2 -= cache[4][2*read_count+1];
			} else {
				lh1 += lgamma(alphas_[ref_allele][i][4][0]+read_count)
						- alphas_[ref_allele][i][4][1];
				lh2 += lgamma(alphas_[ref_allele][i][4][2]+read_count)
						- alphas_[ref_allele][i][4][3];
			}
			log_ret[i] = (lh2 < lh1) ? lh1+log1p(exp(lh2-lh1)) :
			                           lh2+log1p(exp(lh1-lh2)) ;
		}
		return (log_ret - log_ret.maxCoeff()).exp();
	}
	
	DirichletMultinomialMixture(params_t model_a, params_t model_b);
	
protected:
	// NOTE: cache_[a][b][c][d] = sum alpha[a][b][c]+x for x in [0,d/2)
	// NOTE: alpha_[a][b][c] = {alpha1, lgamma(alpha1), alpha2, lgamma(alpha2)}
	typedef std::array<double,2*kCacheSize> cache_type;
	std::array<std::array<std::array<cache_type,5>,10>,5> cache_;
	std::array<std::array<std::array<std::array<double,4>,5>,10>,5> alphas_;
	double f1_, f2_; // log(f) and log(1-f)
};

/*
double DirichletMultinomialLogProbability(double alphas[4], ReadData data) {
	// TODO: Cache most of the math here
	// TODO: Does not include the multinomail coefficient
	int read_count = data.reads[0]+data.reads[1]+data.reads[2]+data.reads[3];
	double alpha_total = alphas[0]+alphas[1]+alphas[2]+alphas[3];
	double result = 0.0;
	for(int i : {0,1,2,3}) {
		for(int x = 0; x < data.reads[i]; ++x) {
			result += log(alphas[i]+x);
		}
	}
	for(int x = 0; x < read_count; ++x)
		result -= log(alpha_total+x);
	return result;
}

DiploidProbs DiploidSequencing(const TetMAParams &params, int ref_allele, ReadData data) {
	DiploidProbs result;
	double alphas_total = (1.0-params.phi_diploid)/params.phi_diploid;
	for(int i : {0,1,2,3}) {
		for(int j=0;j<i;++j) {
			double alphas[4];
			for(int k : {0,1,2,3}) {
				if(k == i || k == j)
					alphas[k] = (0.5-params.error_prob/3.0)*alphas_total;
				else
					alphas[k] = (params.error_prob/3.0)*alphas_total;
			}
			result[i*4+j] = DirichletMultinomialLogProbability(alphas, data);
			result[j*4+i] = result[i*4+j];
		}
		double alphas[4];
		for(int k : {0,1,2,3}) {
			if(k == i)
				alphas[k] = (1.0-params.error_prob)*alphas_total;
			else
				alphas[k] = params.error_prob/3.0*alphas_total;
		}
		result[i*4+i] = DirichletMultinomialLogProbability(alphas, data);
	}
	double scale = result.maxCoeff();
	return (result - scale).exp();
}
*/

} // namespace genotype
} // namespace dng

#endif // DNG_LIKELIHOOD_H

