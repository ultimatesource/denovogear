/*
 * Copyright (c) 2014-2018 Reed A. Cartwright
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
namespace mutation {
/*
  C H I L D
P
A
R  MATRIX
E
N
T
*/

using Matrix = TransitionMatrix;

/* A k-alleles model. See Lewis 2001 and Tuffley and Steel 1997.
   k measures the effective number of alleles, i.e. k-1 is the
   effective number of mutant alleles that can occur at a site.
   This can be calculated from entropy of the mutant-allele spectrum.
   k=5.0 seems to offer a good fit to human segregating variation taking
   into account SNP spectra and indels. 
 */

class Model {
public:
    Model(double u, double k) : u_{u}, k_{k} {
        assert(u_ >= 0.0);
        assert(k_ >= 2.0); 
    }

    Matrix TransitionMatrix(int n);
    Matrix EventTransitionMatrix(int n, int x);
    Matrix MeanTransitionMatrix(int n);

protected:
    double u_;
    double k_;
};


// ret(j,i) = P(i|j)
inline
Matrix Model::TransitionMatrix(int n) {
    assert(n > 0);

    Matrix ret{n,n};
    double beta = u_*k_/(k_-1.0);
    double p_ji = -1.0/k_*expm1(-beta);
    double p_jj = exp(-beta) + p_ji;

    for(int i=0;i<n;++i) {
        for(int j=0;j<n;++j) {
            ret(j,i) = (i == j) ? p_jj : p_ji;
        }
    }
    return ret;    
}

// ret(j,i) = P(i & x mutations | j)
inline
Matrix Model::EventTransitionMatrix(int n, int x) {
    assert(n > 0);
    assert(x >= 0);
    Matrix ret{n,n};
    
    double p_x;
    if(u_ == 0.0) {
        p_x = (x==0) ? 1.0 : 0.0;
    } else {
        p_x = exp(-u_+x*log(u_)-lgamma(x+1));
    }

    double p_ji = (1.0-pow(-1.0/(k_-1.0),x))/k_;
    double p_jj = (1.0+(k_-1.0)*pow(-1.0/(k_-1.0),x))/k_;

    for(int i=0;i<n;++i) {
        for(int j=0;j<n;++j) {
            ret(j,i) = (i == j) ? p_x*p_jj : p_x*p_ji;
        }
    }
    return ret;
}

// ret(j,i) = E[num of mutations | i,j]*P(i|j)
inline
Matrix Model::MeanTransitionMatrix(int n) {
    assert(n > 0);

    Matrix ret{n,n};

    double beta = k_*u_/(k_-1.0);
    double p_jj = -u_/k_*expm1(-beta);
    double p_ji = (u_-p_jj)/(k_-1.0);

    for(int i=0;i<n;++i) {
        for(int j=0;j<n;++j) {
            ret(j,i) = (i == j) ? p_jj : p_ji;
        }
    }
    return ret;
}

struct transition_t {};
struct mean_t {};

inline
Matrix mitosis_haploid_matrix(int size, Model m, transition_t) {
    return m.TransitionMatrix(size);
}

inline
Matrix mitosis_haploid_matrix(int size, Model m, mean_t) {
    return m.MeanTransitionMatrix(size);
}

inline
Matrix mitosis_haploid_matrix(int size, Model m, int count) {
    return m.EventTransitionMatrix(size, count);
}

namespace detail {
inline
void mitosis_diploid_matrix_op(const Matrix& matA, const Matrix& matB, Matrix *p) {
    assert(matA.cols() == matA.rows() && matA.cols() > 0);
    assert(matB.cols() == matB.rows() && matB.cols() > 0);
    assert(matA.cols() == matB.cols());

    const int num_alleles = matB.cols();
    const int num_genotypes = num_alleles*(num_alleles+1)/2;
    assert(p != nullptr && p->cols() == p->rows() && p->cols() == num_genotypes);

    for(int a=0,i=0; a< num_alleles; ++a) { // loop over all child genotypes 
        for(int b=0; b<=a; ++b, ++i) { 
            for(int x=0,j=0; x<num_alleles; ++x) { // loop over all parent genotypes
                for(int y=0; y<=x; ++y,++j) {
                    // x/y => a/b
                    // combine the results from the two phases unless a == b
                    (*p)(j,i) += matA(x,a) * matB(y,b);
                    if(a != b) {
                        (*p)(j,i) += matA(x,b) * matB(y,a);
                    }
                }
            }
        }
    }
}

} // namespace detail

inline
Matrix mitosis_diploid_matrix(int size, Model m, transition_t) {
    assert(size > 0);
    const int num_genotypes = size*(size+1)/2;
    
    Matrix ret = Matrix::Zero(num_genotypes, num_genotypes);

    auto mat = mitosis_haploid_matrix(size, m, transition_t{});
    detail::mitosis_diploid_matrix_op(mat,mat,&ret);

    return ret;
}

inline
Matrix mitosis_diploid_matrix(int size, Model m, int count) {
    assert(size > 0);
    const int num_genotypes = size*(size+1)/2;

    Matrix ret = Matrix::Zero(num_genotypes, num_genotypes);
    
    for(int n=0; n<=count; ++n) {
        auto mat1 = mitosis_haploid_matrix(size, m, n);
        auto mat2 = mitosis_haploid_matrix(size, m, count-n);
        detail::mitosis_diploid_matrix_op(mat1,mat2,&ret);
    }

    return ret;
}

inline
Matrix mitosis_diploid_matrix(int size, Model m, mean_t) {
    assert(size > 0);
    const int num_alleles = size;
    const int num_genotypes = num_alleles*(num_alleles+1)/2;
    
    Matrix ret = Matrix::Zero(num_genotypes, num_genotypes);

    auto mat = mitosis_haploid_matrix(size, m, transition_t{});
    auto avg = mitosis_haploid_matrix(size, m, mean_t{});

    detail::mitosis_diploid_matrix_op(mat, avg, &ret);
    detail::mitosis_diploid_matrix_op(avg, mat, &ret);

    return ret;
}

namespace detail {
inline
void meiosis_haploid_matrix_op(const Matrix& matA, Matrix *p) {
    assert(matA.cols() == matA.rows() && matA.cols() > 0);

    const int num_alleles = matA.cols();
    const int num_genotypes = num_alleles*(num_alleles+1)/2;
    assert(p != nullptr && p->cols() == num_alleles && p->rows() == num_genotypes);

    for(int i=0; i<num_alleles; ++i) { // loop over all child haplotypes
        for(int x=0,j=0; x<num_alleles; ++x) { // loop over all parent genotypes
            for(int y=0; y<=x; ++y,++j) {
                (*p)(j,i) += 0.5*(matA(x,i)+matA(y,i));
            }
        }
    }
}
} // namespace detail

template<typename T>
inline
Matrix meiosis_haploid_matrix(int size, Model m, T arg) {
    assert(size > 0);
    const int num_genotypes = size*(size+1)/2;
    
    Matrix ret = Matrix::Zero(num_genotypes, size);

    auto mat = mitosis_haploid_matrix(size, m, arg);
    detail::meiosis_haploid_matrix_op(mat, &ret);

    return ret;
}

template<typename T>
inline
Matrix mitosis_matrix(int size, Model m, T arg, int parent_ploidy) {
    assert(parent_ploidy == 1 || parent_ploidy == 2);
    if(parent_ploidy == 1) {
        return mitosis_haploid_matrix(size, m, arg);
    }
    return mitosis_diploid_matrix(size, m, arg);    
}

template<typename T>
inline
Matrix gamete_matrix(int size, Model m, T arg, int parent_ploidy) {
    assert(parent_ploidy == 1 || parent_ploidy == 2);
    if(parent_ploidy == 1) {
        return mitosis_haploid_matrix(size, m, arg);
    }
    return meiosis_haploid_matrix(size, m, arg);
}

inline
int number_of_parent_genotypes(const int num_alleles, const int ploidy) {
    assert(ploidy == 1 || ploidy == 2);
    return ((ploidy == 1) ? num_alleles : num_alleles*(num_alleles+1)/2);
}

inline
int number_of_parent_genotype_pairs(const int num_alleles, const int dad_ploidy, const int mom_ploidy) {
    return number_of_parent_genotypes(num_alleles, dad_ploidy)
        * number_of_parent_genotypes(num_alleles, mom_ploidy);
}

namespace detail {
inline
void meiosis_matrix_op(const Matrix& matA, const Matrix& matB, Matrix *p) {
    assert(matA.cols() == matB.cols());

    const int num_alleles = matA.cols();
    const int num_genotypes = num_alleles*(num_alleles+1)/2;
    assert(p != nullptr);
    assert(p->cols() == num_genotypes);
    assert(p->rows() == matA.rows()*matB.rows() );

    for(int a=0,i=0; a < matA.cols(); ++a) { // loop over all child genotypes 
        for(int b=0; b <= a; ++b, ++i) {
            for(int x=0,j=0; x < matA.rows(); ++x) { // loop over all parent genotypes
                for(int y=0; y < matB.rows(); ++y,++j) {
                    // fold the matrix
                    (*p)(j,i) += matA(x,a) * matB(y,b);
                    if(b != a) {
                        (*p)(j,i) += matA(x,b) * matB(y,a);
                    }
                }
            }
        }
    }
}
} // detail

inline
Matrix meiosis_matrix(int size, Model dad_m, Model mom_m, transition_t, int dad_ploidy, int mom_ploidy) {
    assert(dad_ploidy == 1 || dad_ploidy == 2);
    assert(mom_ploidy == 1 || mom_ploidy == 2);

    const int num_alleles = size;
    const int num_genotypes = num_alleles*(num_alleles+1)/2;

    // Construct Mutation Process
    Matrix ret = Matrix::Zero(number_of_parent_genotype_pairs(num_alleles,
            dad_ploidy, mom_ploidy), num_genotypes);
    auto dad = gamete_matrix(size, dad_m, transition_t{}, dad_ploidy);
    auto mom = gamete_matrix(size, mom_m, transition_t{}, mom_ploidy);

    detail::meiosis_matrix_op(dad,mom,&ret);
    return ret;
}

inline
Matrix meiosis_matrix(int size, Model dad_m, Model mom_m, mean_t, int dad_ploidy, int mom_ploidy) {
    assert(dad_ploidy == 1 || dad_ploidy == 2);
    assert(mom_ploidy == 1 || mom_ploidy == 2);

    const int num_alleles = size;
    const int num_genotypes = num_alleles*(num_alleles+1)/2;

    // Construct Mutation Process
    Matrix ret = Matrix::Zero(number_of_parent_genotype_pairs(num_alleles,
            dad_ploidy, mom_ploidy), num_genotypes);

    auto dad = gamete_matrix(size, dad_m, transition_t{}, dad_ploidy);
    auto dad_mean = gamete_matrix(size, dad_m, mean_t{}, dad_ploidy);
    auto mom = gamete_matrix(size, mom_m, transition_t{}, mom_ploidy);
    auto mom_mean = gamete_matrix(size, mom_m, mean_t{}, mom_ploidy);

    detail::meiosis_matrix_op(dad,mom_mean,&ret);
    detail::meiosis_matrix_op(dad_mean,mom,&ret);    
    return ret;
}

inline
Matrix meiosis_matrix(int size, Model dad_m, Model mom_m, int count, int dad_ploidy, int mom_ploidy) {
    assert(dad_ploidy == 1 || dad_ploidy == 2);
    assert(mom_ploidy == 1 || mom_ploidy == 2);

    const int num_alleles = size;
    const int num_genotypes = num_alleles*(num_alleles+1)/2;

    // Construct Mutation Process
    Matrix ret = Matrix::Zero(number_of_parent_genotype_pairs(num_alleles,
            dad_ploidy, mom_ploidy), num_genotypes);

   for(int n=0; n<=count; ++n) {
        auto dad = gamete_matrix(size, dad_m, n, dad_ploidy);
        auto mom = gamete_matrix(size, mom_m, count-n, mom_ploidy);
     
        detail::meiosis_matrix_op(dad,mom,&ret);
    }
    return ret;
}

// k-alleles model from Watterson and Guess (1977) https://doi.org/10.1016/0040-5809(77)90023-5

inline
dng::GenotypeArray population_prior_diploid(int num_obs_alleles, double theta, double hom_bias, double het_bias,
    double kalleles) {
    assert(num_obs_alleles >= 0);

    double k = kalleles;
    double e = theta/(k-1.0);

    double p_hom = (1.0+e)/(1.0+k*e);
    double p_hetk = e/(1.0+k*e);

    double p_RR = p_hom*(2.0+e+(k-1.0)*e*hom_bias)/(2.0+k*e);
    double p_AA = p_hom*(e-e*hom_bias)/(2.0+k*e);

    double p_RA = p_hetk*(2.0+2.0*e+(k-2.0)*e*het_bias)/(2.0+k*e);
    double p_AB = p_hetk*(2.0*e-2.0*e*het_bias)/(2.0+k*e);

    dng::GenotypeArray ret{num_obs_alleles*(num_obs_alleles+1)/2};

    int n=0;
    for(int i=0;i<num_obs_alleles;++i) {
        for(int j=0;j<i;++j) {
            ret(n++) = (j==0 || i==0) ? p_RA : p_AB;
        }
        ret(n++) = (i==0) ? p_RR : p_AA;
    }

    return ret;
}

inline
dng::GenotypeArray population_prior_haploid(int num_obs_alleles, double theta, double hap_bias,
    double kalleles) {
    assert(num_obs_alleles >= 1);

    double k = kalleles;
    double e = theta/(k-1.0);

    double p_R = (1.0+e+(k-1.0)*e*hap_bias)/(1.0+k*e);
    double p_A = (e-e*hap_bias)/(1.0+k*e);

    dng::GenotypeArray ret{num_obs_alleles};
    ret(0) = p_R;
    for(int n=1;n<num_obs_alleles;++n) {
        ret(n) = p_A;
    }

    return ret;
}

inline bool population_prior_check(double theta, double hom_bias, double het_bias, double hap_bias, double kalleles) {
    double k = kalleles;
    double e = theta/(k-1.0);

    if(e < 0) {
        return false;
    }
    if(k < 2.0) {
        return false;
    }
    if(hom_bias > 1.0 || hom_bias < -(2.0+e)/((k-1.0)*e)) {
        return false;
    }
    if(het_bias > 1.0 || het_bias < -(2.0+2.0*e)/((k-2.0)*e)) {
        return false;
    }
    if(hap_bias > 1.0 || hap_bias < -(1.0+e)/((k-1.0)*e)) {
        return false;
    }
    return true;
}

} // namespace mutation

} // namespace dng

#endif // DNG_MUTATION_H
