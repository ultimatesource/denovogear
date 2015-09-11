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
#include <array>

#include <dng/matrix.h>

namespace dng {

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

inline TransitionMatrix mitosis_haploid_matrix(const MutationMatrix &m,
        int mutype = -1) {
    TransitionMatrix ret{4, 4};
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            ret(i, j) = (mutype < 0 || mutype == mitotic_haploid_mutation_counts[i][j]) ?
                        m(i, j) : 0.0;
        }
    }
    return ret; // 4 x 4
}

inline TransitionMatrix mitosis_diploid_matrix(const MutationMatrix &m,
        int mutype = -1) {
    TransitionMatrix ret = TransitionMatrix::Zero(10, 10);
    for(int i = 0; i < 10; ++i) {
        int h = unfolded_diploid_genotypes_upper[i];
        for(int j = 0; j < 16; ++j) {
            int k = folded_diploid_genotypes[j];
            ret(i, k) += (mutype < 0 || mutype == mitotic_diploid_mutation_counts[h][j]) ?
                         kronecker_product(m, m, h, j) : 0.0;
        }
    }
    return ret; // 10 x 10
}

inline TransitionMatrix meiosis_haploid_matrix(const MutationMatrix &m,
        int mutype = -1) {
    TransitionMatrix ret{10, 4};
    auto mm = mitosis_haploid_matrix(m, mutype);
    for(int i = 0; i < 10; ++i) {
        int a = folded_diploid_nucleotides[i][0];
        int b = folded_diploid_nucleotides[i][1];
        for(int j = 0; j < 4; ++j) {
            ret(i, j) = 0.5 * (mm(a, j) + mm(b, j));
        }
    }
    return ret; // 10 x 4
}

inline TransitionMatrix meiosis_diploid_matrix(const MutationMatrix &mdad,
        const MutationMatrix &mmom, int mutype = -1) {
    // Construct Mutation Process
    TransitionMatrix ret = TransitionMatrix::Zero(100, 10);
    TransitionMatrix temp;

    if(mutype <= 0) {
        auto dad = meiosis_haploid_matrix(mdad, mutype);
        auto mom = meiosis_haploid_matrix(mmom, mutype);
        temp = kroneckerProduct(dad, mom);
    } else {
        temp = dng::TransitionMatrix::Zero(100, 16);
        for(int i = 0; i <= mutype; ++i) {
            auto dad = meiosis_haploid_matrix(mdad, i);
            auto mom = meiosis_haploid_matrix(mmom, mutype - i);
            temp += kroneckerProduct(dad, mom);
        }
    }

    // Fold the rows
    for(int i = 0; i < 100; ++i) {
        for(int j = 0; j < 16; ++j) {
            int k = folded_diploid_genotypes[j];
            ret(i, k) += temp(i, j);
        }
    }

    return ret; // 100 x 10
}

inline TransitionMatrix mitosis_haploid_mean_matrix(const MutationMatrix &m) {
    TransitionMatrix ret{4, 4};
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            ret(i, j) = m(i, j) * mitotic_haploid_mutation_counts[i][j];
        }
    }
    return ret; // 4 x 4
}

inline TransitionMatrix mitosis_diploid_mean_matrix(const MutationMatrix &m) {
    TransitionMatrix ret = TransitionMatrix::Zero(10, 10);
    for(int i = 0; i < 10; ++i) {
        int h = unfolded_diploid_genotypes_upper[i];
        for(int j = 0; j < 16; ++j) {
            int k = folded_diploid_genotypes[j];
            ret(i, k) += kronecker_product(m, m, h,
                                           j) * mitotic_diploid_mutation_counts[h][j];
        }
    }
    return ret; // 10 x 10
}

inline TransitionMatrix meiosis_haploid_mean_matrix(const MutationMatrix &m) {
    TransitionMatrix ret{10, 4};
    auto mm = mitosis_haploid_mean_matrix(m);
    for(int i = 0; i < 10; ++i) {
        int a = folded_diploid_nucleotides[i][0];
        int b = folded_diploid_nucleotides[i][1];
        for(int j = 0; j < 4; ++j) {
            ret(i, j) = 0.5 * (mm(a, j) + mm(b, j));
        }
    }
    return ret; // 10 x 4
}

inline TransitionMatrix meiosis_diploid_mean_matrix(const MutationMatrix &mdad,
        const MutationMatrix &mmom) {
    // Construct Mutation Process
    dng::TransitionMatrix ret = dng::TransitionMatrix::Zero(100, 10);

    auto dad = meiosis_haploid_matrix(mdad);
    auto mom = meiosis_haploid_matrix(mmom);
    auto dad_mean = meiosis_haploid_mean_matrix(mdad);
    auto mom_mean = meiosis_haploid_mean_matrix(mmom);

    // Fold the rows
    for(int i = 0; i < 100; ++i) {
        for(int j = 0; j < 16; ++j) {
            int k = folded_diploid_genotypes[j];
            ret(i, k) += kronecker_product(dad, mom_mean, i, j);
            ret(i, k) += kronecker_product(dad_mean, mom, i, j);
        }
    }
    return ret;
}

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


inline dng::GenotypeArray population_prior(double theta,
        const std::array<double, 4> &nuc_freq,
        const std::array<double, 4> &prior) {
    double alpha[4] = {
        theta *nuc_freq[0] + prior[0], theta *nuc_freq[1] + prior[1],
        theta *nuc_freq[2] + prior[2], theta *nuc_freq[3] + prior[3]
    };
    double alpha_sum = alpha[0] + alpha[1] + alpha[2] + alpha[3];
    dng::GenotypeArray ret{10};
    ret <<
        alpha[0]*(1.0 + alpha[0]) / alpha_sum / (1.0 + alpha_sum), // AA
              2.0 * alpha[0]*(alpha[1]) / alpha_sum / (1.0 + alpha_sum), // AC
              2.0 * alpha[0]*(alpha[2]) / alpha_sum / (1.0 + alpha_sum), // AG
              2.0 * alpha[0]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum), // AT
              alpha[1]*(1.0 + alpha[1]) / alpha_sum / (1.0 + alpha_sum), // CC
              2.0 * alpha[1]*(alpha[2]) / alpha_sum / (1.0 + alpha_sum), // CG
              2.0 * alpha[1]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum), // CT
              alpha[2]*(1.0 + alpha[2]) / alpha_sum / (1.0 + alpha_sum), // GG
              2.0 * alpha[2]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum), // GT
              alpha[3]*(1.0 + alpha[3]) / alpha_sum / (1.0 + alpha_sum); // GG
    return ret;
}

} // namespace dng

#endif // DNG_MUTATION_H
