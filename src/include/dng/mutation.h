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

/*
  C H I L D
P
A
R  MATRIX
E
N
T
*/

using MutationMatrix = TransitionMatrix;

namespace Mk {
    // The Mk model of Lewis 2001 and Tuffley and Steel 1997
    // beta = k*a*t and u = (k-1)*a*t
    // num_mutants measures the effective number of alternative alleles (entropy-based)
    //
    inline
    MutationMatrix matrix(int n, double u, double num_mutants) {
        assert(n > 0);
        assert(u >= 0.0);
        assert(num_mutants > 0.0);
        MutationMatrix ret{n,n};
        double K = num_mutants;
        double beta = (K+1.0)*u/K;
        double p_ji = -1.0/(K+1.0)*expm1(-beta);
        double p_jj = exp(-beta) + p_ji;

        for(int i=0;i<n;++i) {
            for(int j=0;j<n;++j) {
                ret(j,i) = (i == j) ? p_jj : p_ji;
            }
        }
        return ret;
    }
} // namespace kalleles

constexpr int MUTATIONS_ALL = -1;
constexpr int MUTATIONS_MEAN = -2;

inline TransitionMatrix mitosis_haploid_matrix(const MutationMatrix &m, const int mutype = MUTATIONS_ALL) {
    assert(m.cols() == m.rows() && m.cols() > 0);

    const int num_alleles = m.cols();

    TransitionMatrix ret{num_alleles, num_alleles};

    for(int i = 0; i < num_alleles; ++i) { // loop over children
        for(int j = 0; j < num_alleles; ++j) { // loop over parents
            int u = (i != j);
            if(mutype != MUTATIONS_MEAN) {
                // Update as needed
                u = (mutype == MUTATIONS_ALL || u == mutype) ? 1 : 0;
            }
            ret(j,i) = m(j,i)*u;
        }
    }
    return ret;
}

inline TransitionMatrix mitosis_diploid_matrix(const MutationMatrix &m, const int mutype = MUTATIONS_ALL) {
    assert(m.cols() == m.rows() && m.cols() > 0);
    
    const int num_alleles = m.cols();
    const int num_genotypes = num_alleles*(num_alleles+1)/2;

    TransitionMatrix ret{num_genotypes,num_genotypes};

    for(int a=0,i=0; a<num_alleles; ++a) { // loop over all child genotypes 
        for(int b=0; b<=a; ++b, ++i) { 
            for(int x=0,j=0; x<num_alleles; ++x) { // loop over all parent genotypes
                for(int y=0; y<=x; ++y,++j) {
                    // x/y => a/b
                    int u = (x != a) + (y != b); // number of mutations for phase 1
                    int v = (x != b) + (y != a); // number of mutations for phase 2
                    if(mutype != MUTATIONS_MEAN) {
                        // Update u and v as needed
                        u = (mutype == MUTATIONS_ALL || u == mutype) ? 1 : 0;
                        v = (mutype == MUTATIONS_ALL || v == mutype) ? 1 : 0;
                    }
                    // combine the results from the two phases unless a == b
                    ret(j,i) = m(x,a)*m(y,b)*u + m(x,b)*m(y,a)*((a!=b) ? v : 0);
                }
            }
        }
    }
    return ret;
}

inline TransitionMatrix meiosis_haploid_matrix(const MutationMatrix &m, int mutype = MUTATIONS_ALL) {
    assert(m.cols() == m.rows() && m.cols() > 0);
    
    const int num_alleles = m.cols();
    const int num_genotypes = num_alleles*(num_alleles+1)/2;

    TransitionMatrix ret{num_genotypes, num_alleles};

    for(int i=0; i<num_alleles; ++i) { // loop over all child haplotypes
        for(int x=0,j=0; x<num_alleles; ++x) { // loop over all parent genotypes
            for(int y=0; y<=x; ++y,++j) {
                // x/y => i
                int u = (i != x); // number of mutations from chrom 1
                int v = (i != y); // number of mutations from chrom 2
                if(mutype != MUTATIONS_MEAN) {
                    // Update u and v as needed
                    u = (mutype == MUTATIONS_ALL || u == mutype) ? 1 : 0;
                    v = (mutype == MUTATIONS_ALL || v == mutype) ? 1 : 0;
                }
                ret(j,i) = 0.5*(m(x,i)*u + m(y,i)*v);
            }
        }
    }
    return ret;
}

inline TransitionMatrix mitosis_matrix(const int parent_ploidy, const MutationMatrix &m, const int mutype = MUTATIONS_ALL) {
    assert(parent_ploidy == 1 || parent_ploidy == 2);
    return (parent_ploidy == 1) ? mitosis_haploid_matrix(m, mutype) : mitosis_diploid_matrix(m,mutype);    
}

inline TransitionMatrix gamete_matrix(const int parent_ploidy, const MutationMatrix &m, const int mutype = MUTATIONS_ALL) {
    assert(parent_ploidy == 1 || parent_ploidy == 2);
    return (parent_ploidy == 1) ? mitosis_haploid_matrix(m, mutype) : meiosis_haploid_matrix(m,mutype);
}

inline int number_of_parent_genotypes(const int num_alleles, const int ploidy) {
    assert(ploidy == 1 || ploidy == 2);
    return ((ploidy == 1) ? num_alleles : num_alleles*(num_alleles+1)/2);
}

inline int number_of_parent_genotype_pairs(const int num_alleles, const int dad_ploidy, const int mom_ploidy) {
    return number_of_parent_genotypes(num_alleles, dad_ploidy)
        * number_of_parent_genotypes(num_alleles, mom_ploidy);
}

inline TransitionMatrix meiosis_matrix(const int dad_ploidy, const MutationMatrix &dad_m,
    const int mom_ploidy, const MutationMatrix &mom_m, const int mutype = MUTATIONS_ALL) {
    assert(dad_ploidy == 1 || dad_ploidy == 2);
    assert(mom_ploidy == 1 || mom_ploidy == 2);
    assert(dad_m.cols() == dad_m.rows() && dad_m.cols() > 0);
    assert(mom_m.cols() == mom_m.rows() && mom_m.cols() > 0);
    assert(mom_m.cols() == dad_m.cols());
    assert(mutype >= 0 || mutype == MUTATIONS_ALL || mutype == MUTATIONS_MEAN);

    const int num_alleles = dad_m.cols();
    const int num_genotypes = num_alleles*(num_alleles+1)/2;

    // Construct Mutation Process
    TransitionMatrix temp{number_of_parent_genotype_pairs(num_alleles,
            dad_ploidy, mom_ploidy), num_alleles*num_alleles};
    if(mutype == MUTATIONS_ALL) {
        TransitionMatrix dad = gamete_matrix(dad_ploidy, dad_m, MUTATIONS_ALL);
        TransitionMatrix mom = gamete_matrix(mom_ploidy, mom_m, MUTATIONS_ALL);
        temp = kroneckerProduct(dad, mom);
    } else if(mutype == MUTATIONS_MEAN) {
        TransitionMatrix dad = gamete_matrix(dad_ploidy, dad_m, MUTATIONS_ALL);
        TransitionMatrix dad_mean = gamete_matrix(dad_ploidy, dad_m, MUTATIONS_MEAN);
        TransitionMatrix mom = gamete_matrix(mom_ploidy, mom_m, MUTATIONS_ALL);
        TransitionMatrix mom_mean = gamete_matrix(mom_ploidy, mom_m, MUTATIONS_MEAN);
        temp = kroneckerProduct(dad,mom_mean)+kroneckerProduct(dad_mean,mom);
    } else {
        temp.setZero();
        for(int i = 0; i <= mutype; ++i) {
            TransitionMatrix dad = gamete_matrix(dad_ploidy, dad_m, i);
            TransitionMatrix mom = gamete_matrix(mom_ploidy, mom_m, mutype - i);
            temp += kroneckerProduct(dad, mom);
        }
    }

    // Fold the rows
    TransitionMatrix ret{temp.rows(), num_genotypes};
    for(int a=0,i=0; a<num_alleles; ++a) { // loop over all child genotypes 
        for(int b=0; b<=a; ++b, ++i) {
            for(int j=0; j<temp.rows(); ++j) {
                // fold the temp matrix
                ret(j,i) = temp(j, a*num_alleles+b);
                if(b != a) {
                    ret(j,i) += temp(j, b*num_alleles+a);
                }
            }
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
    for(int n=1;n<=num_alts;++n) {
        ret(n) = p_A;
    }

    return ret;
}

} // namespace dng

#endif // DNG_MUTATION_H
