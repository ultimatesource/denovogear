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
#ifndef DNG_MUTATION_H
#define DNG_MUTATION_H


#include <cmath>
#include <cassert>
#include <array>
#include <iostream>
#include <dng/matrix.h>
#include <dng/genotyper.h>

namespace dng {

typedef Eigen::Matrix4d MutationMatrix;

namespace f81 {
// Construct an F81 mutation matrix
inline MutationMatrix matrix(double mu, const std::array<double, 4> &nuc_freq) {
    double beta = 1.0;
    for(auto d : nuc_freq) {
        beta -= d * d;
    }
    double p = -expm1(-mu / beta);

    MutationMatrix ret;
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            ret(i, j) = nuc_freq[j] * p;
        }
        ret(i, i) += 1.0 - p;
    }
    return ret;
}

// Calculate the maximum likelihood estimate of mu given q fraction of observed differences
inline double estimate(double q, const std::array<double, 4> &nuc_freq) {
    double beta = 1.0;
    for(auto d : nuc_freq) {
        beta -= d * d;
    }
    return -beta * log(1.0 - q / beta);
}

} // namespace f81

namespace Mk {
    // The Mk model of Lewis 2001 and Tuffley and Steel 1997
    // u = (k-1)*a*t

    MutationMatrix matrix(double u, int k) {
        MutationMatrix ret{k,k};
        double K = k;
        double beta = K*u/(K-1.0);
        double p_ij = -1.0/K*expm1(-beta);
        double p_ii = exp(-beta) + p_ij;

        for(int i=0;i<k;++i) {
            for(int j=0;j<k;++j) {
                ret(j,i) = (i == j) ? p_ii : p_ij;
            }
        }
        return ret;
    }

} // namespace kalleles

constexpr int MUTATIONS_ALL = -1;
constexpr int MUTATIONS_MEAN = -2;

constexpr int mitotic_haploid_mutation_counts[4][4] = {
    {0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}
};

constexpr char mitotic_haploid_mutation_labels[4][4][4] = {
    {"A>A", "A>C", "A>G", "A>T"}, {"C>A", "C>C", "C>G", "C>T"},
    {"G>A", "G>C", "G>G", "G>T"}, {"T>A", "T>C", "T>G", "T>T"},
};

constexpr int mitotic_diploid_mutation_counts[16][16] = {
    {0, 1, 1, 1, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2},
    {1, 0, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2},
    {1, 1, 0, 1, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2},
    {1, 1, 1, 0, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1},
    {1, 2, 2, 2, 0, 1, 1, 1, 1, 2, 2, 2, 1, 2, 2, 2},
    {2, 1, 2, 2, 1, 0, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2},
    {2, 2, 1, 2, 1, 1, 0, 1, 2, 2, 1, 2, 2, 2, 1, 2},
    {2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 1, 2, 2, 2, 1},
    {1, 2, 2, 2, 1, 2, 2, 2, 0, 1, 1, 1, 1, 2, 2, 2},
    {2, 1, 2, 2, 2, 1, 2, 2, 1, 0, 1, 1, 2, 1, 2, 2},
    {2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 0, 1, 2, 2, 1, 2},
    {2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 1},
    {1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 0, 1, 1, 1},
    {2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 1, 0, 1, 1},
    {2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 0, 1},
    {2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 0}
};

extern const char mitotic_diploid_mutation_labels[10][10][6];
extern const char meiotic_diploid_mutation_labels[100][10][9];

inline TransitionMatrix mitosis_haploid_matrix(const MutationMatrix &m, const int mutype = MUTATIONS_ALL) {
    TransitionMatrix ret{4, 4};
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            if(mutype == MUTATIONS_MEAN) {
                ret(i, j) = m(i,j) * mitotic_haploid_mutation_counts[i][j];
            } else if( mutype == MUTATIONS_ALL || mutype == mitotic_haploid_mutation_counts[i][j] ) {
                ret(i,j) = m(i,j);
            } else {
                ret(i,j) = 0.0;
            }
        }
    }
    return ret; // 4 x 4
}

inline TransitionMatrix mitosis_diploid_matrix(const MutationMatrix &m, const int mutype = MUTATIONS_ALL) {
    TransitionMatrix ret = TransitionMatrix::Zero(10, 10);
    for(int i = 0; i < 10; ++i) {
        int h = genotype::unfolded_diploid_genotypes_upper[i];
        for(int j = 0; j < 16; ++j) {
            int k = genotype::folded_diploid_genotypes[j];
            if(mutype == MUTATIONS_MEAN) {
                ret(i, k) += kronecker_product_coef(m, m, h, j) * mitotic_diploid_mutation_counts[h][j];
            } else if( mutype == MUTATIONS_ALL || mutype == mitotic_diploid_mutation_counts[h][j] ) {
                ret(i,k) += kronecker_product_coef(m, m, h, j);
            } else {
                ret(i,k) += 0.0;
            }
        }
    }
    return ret; // 10 x 10
}

inline TransitionMatrix meiosis_haploid_matrix(const MutationMatrix &m, int mutype = MUTATIONS_ALL) {
    TransitionMatrix ret{10, 4};
    auto mm = mitosis_haploid_matrix(m, mutype);
    for(int i = 0; i < 10; ++i) {
        int a = genotype::folded_diploid_nucleotides[i][0];
        int b = genotype::folded_diploid_nucleotides[i][1];
        for(int j = 0; j < 4; ++j) {
            ret(i, j) = 0.5 * (mm(a, j) + mm(b, j));
        }
    }
    return ret; // 10 x 4
}

inline TransitionMatrix mitosis_matrix(const int parent_ploidy, const MutationMatrix &m, const int mutype = MUTATIONS_ALL) {
    assert(parent_ploidy == 1 || parent_ploidy == 2);
    return (parent_ploidy == 1) ? mitosis_haploid_matrix(m, mutype) : mitosis_diploid_matrix(m,mutype);    
}

inline TransitionMatrix gamete_matrix(const int parent_ploidy, const MutationMatrix &m, const int mutype = MUTATIONS_ALL) {
    assert(parent_ploidy == 1 || parent_ploidy == 2);
    return (parent_ploidy == 1) ? mitosis_haploid_matrix(m, mutype) : meiosis_haploid_matrix(m,mutype);
}

inline int number_of_parent_genotypes(const int ploidy) {
    assert(ploidy == 1 || ploidy == 2);
    return ((ploidy == 1) ? 4 : 10);
}

inline int number_of_parent_genotype_pairs(const int dad_ploidy, const int mom_ploidy) {
    return number_of_parent_genotypes(dad_ploidy)*number_of_parent_genotypes(mom_ploidy);
}

inline TransitionMatrix meiosis_matrix(const int dad_ploidy, const MutationMatrix &dad_m,
    const int mom_ploidy, const MutationMatrix &mom_m, const int mutype = MUTATIONS_ALL) {
    assert(dad_ploidy == 1 || dad_ploidy == 2);
    assert(mom_ploidy == 1 || mom_ploidy == 2);

    // Construct Mutation Process
    TransitionMatrix temp;
    if(mutype == MUTATIONS_MEAN) {
        TransitionMatrix dad = gamete_matrix(dad_ploidy, dad_m, MUTATIONS_ALL);
        TransitionMatrix dad_mean = gamete_matrix(dad_ploidy, dad_m, MUTATIONS_MEAN);
        TransitionMatrix mom = gamete_matrix(mom_ploidy, mom_m, MUTATIONS_ALL);
        TransitionMatrix mom_mean = gamete_matrix(mom_ploidy, mom_m, MUTATIONS_MEAN);
        temp = kroneckerProduct(dad,mom_mean)+kroneckerProduct(dad_mean,mom);
    } else if(mutype == MUTATIONS_ALL) {
        TransitionMatrix dad = gamete_matrix(dad_ploidy, dad_m, MUTATIONS_ALL);
        TransitionMatrix mom = gamete_matrix(mom_ploidy, mom_m, MUTATIONS_ALL);
        temp = kroneckerProduct(dad, mom);
    } else if(mutype >= 0) {
        temp = dng::TransitionMatrix::Zero(number_of_parent_genotype_pairs(dad_ploidy,mom_ploidy), 16);
        for(int i = 0; i <= mutype; ++i) {
            TransitionMatrix dad = gamete_matrix(dad_ploidy, dad_m, i);
            TransitionMatrix mom = gamete_matrix(mom_ploidy, mom_m, mutype - i);
            temp += kroneckerProduct(dad, mom);
        }
    } else {
        assert(false); // should not get here
    }

    // Fold the rows
    TransitionMatrix ret = TransitionMatrix::Zero(temp.rows(), 10);
    for(int i = 0; i < temp.rows(); ++i) {
        for(int j = 0; j < temp.cols(); ++j) {
            int k = genotype::folded_diploid_genotypes[j];
            ret(i, k) += temp(i, j);
        }
    }

    return ret;
}

inline TransitionMatrix meiosis_diploid_mean_matrix(const MutationMatrix &mdad,
        const MutationMatrix &mmom) {
    return meiosis_matrix(2,mdad,2,mmom, MUTATIONS_MEAN);
}

inline TransitionMatrix meiosis_diploid_matrix(const MutationMatrix &mdad,
        const MutationMatrix &mmom, int mutype = MUTATIONS_ALL) {    
    return meiosis_matrix(2,mdad,2,mmom,mutype);
}

inline TransitionMatrix meiosis_diploid_matrix_xlinked(const MutationMatrix &mdad,
        const MutationMatrix &mmom, int mutype = MUTATIONS_ALL) {
    return meiosis_matrix(1,mdad,2,mmom,mutype);
}

inline
std::array<double, 4> population_alphas(double theta,
        const std::array<double, 4> &nuc_freq,
        const std::array<double, 4> &prior) {
    assert(nuc_freq[0] >= 0 && nuc_freq[1] >= 0 && nuc_freq[2] >= 0 && nuc_freq[3] >= 0);
    assert(prior[0] >= 0 && prior[1] >= 0 && prior[2] >= 0 && prior[3] >= 0);

    double nuc_sum = nuc_freq[0] + nuc_freq[1] + nuc_freq[2] + nuc_freq[3];
    theta = theta/nuc_sum;

    return {{theta*nuc_freq[0] + prior[0], theta*nuc_freq[1] + prior[1],
             theta*nuc_freq[2] + prior[2], theta*nuc_freq[3] + prior[3]
    }};
}

inline dng::GenotypeArray population_prior_diploid_dm(double theta,
        const std::array<double, 4> &nuc_freq,
        const std::array<double, 4> &prior) {
    std::array<double, 4> alpha = population_alphas(theta, nuc_freq, prior);

    double alpha_sum = alpha[0] + alpha[1] + alpha[2] + alpha[3];
    dng::GenotypeArray ret{10};
    for(int i=0;i<10;++i) {
        int n1 = genotype::folded_diploid_nucleotides[i][0];
        int n2 = genotype::folded_diploid_nucleotides[i][1];
        if(n1 == n2) {
            ret(i) = alpha[n1]*(1.0 + alpha[n1]) / alpha_sum / (1.0 + alpha_sum);
        } else {
            ret(i) = 2.0 * alpha[n1]*(alpha[n2]) / alpha_sum / (1.0 + alpha_sum);

        }
    }
    return ret;
}


inline dng::GenotypeArray population_prior_haploid_dm(double theta,
        const std::array<double, 4> &nuc_freq,
        const std::array<double, 4> &prior) {
    std::array<double, 4> alpha = population_alphas(theta, nuc_freq, prior);

    double alpha_sum = alpha[0] + alpha[1] + alpha[2] + alpha[3];
    dng::GenotypeArray ret{4};
    ret <<  alpha[0] / alpha_sum,
            alpha[1] / alpha_sum,
            alpha[2] / alpha_sum,
            alpha[3] / alpha_sum;
    return ret;
}

inline bool population_prior_check_ia(double theta, double hom_bias, double het_bias, double hap_bias) {
    return (theta >= 0)
        && (theta*hom_bias >= -2.0 && hom_bias <= 1.0)
        && (theta*het_bias >= -2.0 && het_bias <= 1.0)
        && (theta*hap_bias >= -1.0 && hap_bias <= 1.0)
        ;
}

inline
dng::GenotypeArray population_prior_diploid_ia(double theta, double hom_bias, double het_bias,
    int num_alts, bool known_anc=true) {
    assert(num_alts >= 0);

    // Check to see if we have any alleles
    if(!known_anc && num_alts == 0) {
        return {};
    }

    double p_hom = 1.0/(1.0+theta);
    double p_het = theta/(1.0+theta);

    double k = num_alts;

    double p_RR=0.0, p_AA=0.0, p_RA=0.0, p_AB=0.0;
    if(known_anc) {
        if(num_alts == 0) {
            p_RR = 1.0;
        } else {
            p_RR = p_hom*(2.0+theta*hom_bias)/(2.0+theta);
            p_AA = p_hom*theta*(1.0-hom_bias)/(2.0+theta)*(1.0/k);
            if(num_alts == 1) {
                p_RA = p_het;
            } else {
                p_RA = p_het*(2.0+theta*het_bias)/(2.0+theta)*(1.0/k);
                p_AB = p_het*theta*(1.0-het_bias)/(2.0+theta)*(2.0/(k*(k-1.0)));
            }
        }
    } else {
        if(num_alts == 1) {
            p_AA = 1.0;
        } else {
            p_AA = p_hom*(1.0/k);
            p_AB = p_het*(2.0/(k*(k-1.0)));
        }
    }

    dng::GenotypeArray ret{(num_alts+1)*(num_alts+2)/2};

    int n=0;
    for(int i=0;i<(num_alts+1);++i) {
        for(int j=0;j<i;++j) {
            ret(n++) = (j==0 || i==0) ? p_RA : p_AB;
        }
        ret(n++) = (i==0) ? p_RR : p_AA;
    }

    return ret;
}

inline
dng::GenotypeArray population_prior_haploid_ia(double theta, double hap_bias,
    int num_alts, bool known_anc=true) {
    assert(num_alts >= 0);

    // Check to see if we have any alleles
    if(!known_anc && num_alts == 0) {
        return {};
    }

    double k = num_alts;
    double p_R = 0.0, p_A = 0.0;
    if(known_anc) {
        if(num_alts == 0) {
            p_R = 1.0;
        } else {
            p_R = (1.0+theta*hap_bias)/(1.0+theta);
            p_A = theta*(1.0-hap_bias)/(1.0+theta)*(1.0/k);
        }
    } else {
        p_A = 1.0/k;
    }
    
    dng::GenotypeArray ret{num_alts+1};
    ret(0) = p_R;
    for(int n=0;n<(num_alts+1);++n) {
        ret(n) = p_A;
    }

    return ret;
}

} // namespace dng

#endif // DNG_MUTATION_H
