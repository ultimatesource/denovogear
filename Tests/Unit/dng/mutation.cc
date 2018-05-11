/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2017-2018 Reed A. Cartwright
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
 *           Reed A. Cartwright <reed@cartwrig.ht>
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

#define BOOST_TEST_MODULE dng::mutation

#include <dng/mutation.h>

#include <dng/detail/rangeio.h>

#include <iostream>
#include <iomanip>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/multi_array.hpp>

#include "../testing.h"

#include "../boost_test_helper.h"

using namespace dng;
using namespace dng::genotype;
using namespace dng::mutation;

using dng::detail::make_test_range;
using d4 = std::array<double,4>;

// This generic function allows us to test the same parameters
// for different test functions
template<typename F>
void run_mutation_tests(F test, double prec = 2*DBL_EPSILON) {
    test(4, 0.0,  4, prec);
    // test(1, 1e-8, 4, prec);
    
    // test(2, 1e-8, 4, prec);
    // test(3, 1e-8, 4, prec);
    // test(4, 1e-8, 4, prec);
    
    // test(4, 1e-8, 4, prec);
    // test(4, 1e-9, 5, prec);
    // test(4, 1e-6, 6, prec);
    // test(4, 1e-3, 7, prec);
}

BOOST_AUTO_TEST_CASE(test_model) {
    auto test = [](int n, double u, double k, double prec) -> void {
        BOOST_TEST_CONTEXT("n=" << n << ", u=" << u << ", k=" << k)
    {
        assert(n <= k);

        using namespace boost::numeric::ublas;
        using mat = matrix<double,column_major,std::vector<double>>;
        
        // Build matrix
        mat Q(n+1,n+1);

        vector<double> freqs(n+1, 1.0/k);
        freqs[n] = 1.0-n/k;

        for(int i=0; i<=n; ++i) {
            for(int j=0; j<=n; ++j) {
                Q(i,j) = freqs[j];
            }
            Q(i,i) = Q(i,i)-1.0;
        }

        // Scale matrix into substitution time
        Q *= k/(k-1);

        mat P = identity_matrix<double>(n+1);
        mat m = P;
        double t = 1.0;
        double f = 1.0;
        // Use the first 11 elements of the expm series 
        for(int g=1;g <= 11; ++g) {
            mat a = prec_prod(m,Q);
            m = a;
            t *= u;
            f *= g;
            P += m * (t/f);
        }

        std::vector<double> expected_matrix;
        for(int i=0;i<n;++i) {
            for(int j=0;j<n;++j) {
                expected_matrix.push_back(P(j,i));
            }
        }

        auto model = Model{u,k};
        auto test_matrix_ = model.TransitionMatrix(n);
        auto test_matrix = make_test_range(test_matrix_);
        CHECK_CLOSE_RANGES( test_matrix, expected_matrix, prec);
    }
    };

    run_mutation_tests(test);
}


BOOST_AUTO_TEST_CASE(test_mitosis_haploid_matrix) {
    auto test = [](int n, double u, double k, double prec) -> void {
        BOOST_TEST_CONTEXT("n=" << n << ", u=" << u << ", k=" << k)
    {
        auto m = Model{u,k}.TransitionMatrix(n);

        // Test transitions
        auto x_all = mitosis_haploid_matrix(n, Model{u,k}, transition_t{});
        auto expected_all = make_test_range(m);
        auto test_all = make_test_range(x_all);
        CHECK_CLOSE_RANGES(test_all, expected_all, prec);

        // Test 0 Mutations
        auto x_0 = mitosis_haploid_matrix(n, Model{u,k}, 0);
        auto m_0 = m;
        for(int i=0;i<n;++i) {
            for(int j=0;j<n;++j) {
                if(i != j) {
                    m_0(i,j) = 0.0;
                }
            }
        }
        auto expected_0 = make_test_range(m_0);
        auto test_0 = make_test_range(x_0);
        CHECK_CLOSE_RANGES(test_0, expected_0, prec);

        // Test 1 Mutation
        auto x_1 = mitosis_haploid_matrix(n, Model{u,k}, 1);
        auto m_1 = m;
        for(int i=0;i<n;++i) {
            for(int j=0;j<n;++j) {
                if(i == j) {
                    m_1(i,j) = 0.0;
                }
            }
        }
        auto expected_1 = make_test_range(m_1);
        auto test_1 = make_test_range(x_1);
        CHECK_CLOSE_RANGES(test_1, expected_1, prec);

        // Test 2 Mutations
        auto x_2 = mitosis_haploid_matrix(n, Model{u,k}, 2);
        auto m_2 = m;
        for(int i=0;i<n;++i) {
            for(int j=0;j<n;++j) {
                m_2(i,j) = 0.0;
            }
        }
        auto expected_2 = make_test_range(m_2);
        auto test_2 = make_test_range(x_2);
        CHECK_CLOSE_RANGES(test_2, expected_2, prec);

        // Test Mean Mutations
        auto x_mean = mitosis_haploid_matrix(n, Model{u,k}, mean_t{});
        auto m_mean = m;
        for(int i=0;i<n;++i) {
            for(int j=0;j<n;++j) {
                if(i == j) {
                    m_mean(i,j) = 0.0;
                }
            }
        }
        auto expected_mean = make_test_range(m_mean);
        auto test_mean = make_test_range(x_mean);
        CHECK_CLOSE_RANGES(test_mean, expected_mean, prec);
    }};

    run_mutation_tests(test);
}


BOOST_AUTO_TEST_CASE(test_mitosis_diploid_matrix) {
    auto test = [](int n, double u, double k, double prec) -> void {
        BOOST_TEST_CONTEXT("n=" << n << ", u=" << u << ", k=" << k)
    {
        auto f = Model{u,k}.TransitionMatrix(n);
        int n2 = (n*(n+1)/2);

        // compute n*n x n*n mutation matrix
        boost::multi_array<double, 4> mk(boost::extents[n][n][n][n]);

        for(int a=0;a<n;++a) {
            for(int b=0;b<n;++b) {
                for(int x=0;x<n;++x) {
                    for(int y=0;y<n;++y) {
                        // a -> x and b -> y
                        mk[a][b][x][y] = f(a,x)*f(b,y);
                    }
                }
            }
        }
        // fold the matrix
        std::vector<double> m( n2*n2 );
        for(int x=0,z=0;x<n;++x) {
            for(int y=0;y<=x;++y) {
                for(int a=0,c=0;a<n;++a) {
                    for(int b=0;b<=a;++b) {
                        m[z] = mk[a][b][x][y];
                        if( x != y) {
                            m[z] += mk[a][b][y][x];
                        }
                        ++z;
                    }
                }
            }
        }

        // Test MUTATIONS_ALL
        auto x_all = mitosis_diploid_matrix(n, Model{u,k}, transition_t{});
        auto expected_all = m;
        auto test_all = make_test_range(x_all);
        CHECK_CLOSE_RANGES(test_all, expected_all, prec);

        // Test 0 Mutations
        auto x_0 = mitosis_diploid_matrix(n, Model{u,k}, 0);
        for(int x=0,z=0;x<n;++x) {
            for(int y=0;y<=x;++y) {
                for(int a=0;a<n;++a) {
                    for(int b=0;b<=a;++b) {
                        if(a == x && b == y) {
                            m[z] = mk[a][b][x][y];
                        } else {
                            m[z] = 0.0;
                        }
                        ++z;
                    }
                }
            }
        }
        auto expected_0 = m;
        auto test_0 = make_test_range(x_0);
        CHECK_CLOSE_RANGES(test_0, expected_0, prec);

        // Test 1 Mutations
        auto x_1 = mitosis_diploid_matrix(n, Model{u,k}, 1);
        for(int x=0,z=0;x<n;++x) {
            for(int y=0;y<=x;++y) {
                for(int a=0;a<n;++a) {
                    for(int b=0;b<=a;++b) {
                        m[z] = 0.0;
                        if(((a != x) + (b != y)) == 1) {
                            m[z] += mk[a][b][x][y];
                        }
                        if(x != y && ((a != y) + (b != x)) == 1) {
                            m[z] += mk[a][b][y][x];
                        }
                        ++z;
                    }
                }
            }
        }
        auto expected_1 = m;
        auto test_1 = make_test_range(x_1);
        CHECK_CLOSE_RANGES(test_1, expected_1, prec);

        // Test 2 Mutations
        auto x_2 = mitosis_diploid_matrix(n, Model{u,k}, 2);
        for(int x=0,z=0;x<n;++x) {
            for(int y=0;y<=x;++y) {
                for(int a=0;a<n;++a) {
                    for(int b=0;b<=a;++b) {
                        m[z] = 0.0;
                        if(((a != x) + (b != y)) == 2) {
                            m[z] += mk[a][b][x][y];
                        }
                        if(x != y && ((a != y) + (b != x)) == 2) {
                            m[z] += mk[a][b][y][x];
                        }
                        ++z;
                    }
                }
            }
        }
        auto expected_2 = m;
        auto test_2 = make_test_range(x_2);
        CHECK_CLOSE_RANGES(test_2, expected_2, prec);

        // Test 3 Mutations
        auto x_3 = mitosis_diploid_matrix(n, Model{u,k}, 3);
        for(int x=0,z=0;x<n;++x) {
            for(int y=0;y<=x;++y) {
                for(int a=0;a<n;++a) {
                    for(int b=0;b<=a;++b) {
                        m[z] = 0.0;
                        ++z;
                    }
                }
            }
        }
        auto expected_3 = m;
        auto test_3 = make_test_range(x_3);
        CHECK_CLOSE_RANGES(test_3, expected_3, prec);

        // Test Mean Mutations
        auto x_mean = mitosis_diploid_matrix(n, Model{u,k}, mean_t{});
        for(int x=0,z=0;x<n;++x) {
            for(int y=0;y<=x;++y) {
                for(int a=0;a<n;++a) {
                    for(int b=0;b<=a;++b) {
                        m[z] = 0.0;
                        m[z] += mk[a][b][x][y]*((a != x) + (b != y));
                        if(x != y) {
                            m[z] += mk[a][b][y][x]*((a != y) + (b != x));
                        }
                        ++z;
                    }
                }
            }
        }
        auto expected_mean = m;
        auto test_mean = make_test_range(x_mean);
        CHECK_CLOSE_RANGES(test_mean, expected_mean, prec);
    }
    };

    run_mutation_tests(test);
}


BOOST_AUTO_TEST_CASE(test_meiosis_haploid_matrix) {
    auto test = [](int n, double u, double k, double prec) -> void {
        BOOST_TEST_CONTEXT("n=" << n << ", u=" << u << ", k=" << k)
    {
        auto f = Model{u,k}.TransitionMatrix(n);;
        int n2 = (n*(n+1)/2);
        
        std::vector<double> m(n*n2);
        for(int x=0,z=0;x<n;++x) {
            for(int a=0,c=0;a<n;++a) {
                for(int b=0;b<=a;++b) {
                    m[z] = 0.5*(f(a,x)+f(b,x));
                    ++z;
                }
            }
        }

        // Test MUTATIONS_ALL
        auto x_all = meiosis_haploid_matrix(n, Model{u,k}, transition_t{});
        auto expected_all = m;
        auto test_all = make_test_range(x_all);
        CHECK_CLOSE_RANGES(test_all, expected_all, prec);

        // Test 0 Mutations
        auto x_0 = meiosis_haploid_matrix(n, Model{u,k}, 0);
        for(int x=0,z=0;x<n;++x) {
            for(int a=0,c=0;a<n;++a) {
                for(int b=0;b<=a;++b) {
                    m[z] = 0.0;
                    m[z] += 0.5*f(a,x)*(a==x);
                    m[z] += 0.5*f(b,x)*(b==x);
                    ++z;
                }
            }
        }
        auto expected_0 = m;
        auto test_0 = make_test_range(x_0);
        CHECK_CLOSE_RANGES(test_0, expected_0, prec);

        // Test 1 Mutation
        auto x_1 = meiosis_haploid_matrix(n, Model{u,k}, 1);
        for(int x=0,z=0;x<n;++x) {
            for(int a=0,c=0;a<n;++a) {
                for(int b=0;b<=a;++b) {
                    m[z] = 0.0;
                    m[z] += 0.5*f(a,x)*(a!=x);
                    m[z] += 0.5*f(b,x)*(b!=x);
                    ++z;
                }
            }
        }
        auto expected_1 = m;
        auto test_1 = make_test_range(x_1);
        CHECK_CLOSE_RANGES(test_1, expected_1, prec);

        // Test 2 Mutations
        auto x_2 = meiosis_haploid_matrix(n, Model{u,k}, 2);
        for(int x=0,z=0;x<n;++x) {
            for(int a=0,c=0;a<n;++a) {
                for(int b=0;b<=a;++b) {
                    m[z] = 0.0;
                    ++z;
                }
            }
        }
        auto expected_2 = m;
        auto test_2 = make_test_range(x_2);
        CHECK_CLOSE_RANGES(test_2, expected_2, prec);

        // Test Mean Mutations
        auto x_mean = meiosis_haploid_matrix(n, Model{u,k}, mean_t{});
        for(int x=0,z=0;x<n;++x) {
            for(int a=0,c=0;a<n;++a) {
                for(int b=0;b<=a;++b) {
                    m[z] = 0.0;
                    m[z] += 0.5*f(a,x)*(a!=x);
                    m[z] += 0.5*f(b,x)*(b!=x);
                    ++z;
                }
            }
        }
        auto expected_mean = m;
        auto test_mean = make_test_range(x_mean);
        CHECK_CLOSE_RANGES(test_mean, expected_mean, prec);

    }};

    run_mutation_tests(test);
}


BOOST_AUTO_TEST_CASE(test_mitosis_matrix) {
    auto test = [](int n, double u, double k, double prec) -> void {
        BOOST_TEST_CONTEXT("n=" << n << ", u=" << u << ", k=" << k)
    {
        auto f = Model{u,k}.TransitionMatrix(n);;
        int n2 = (n*(n+1)/2);

        // Ploidy 1
        BOOST_TEST_CONTEXT("ploidy=1") {
            BOOST_TEST_CONTEXT("mutype=all") {
                auto m = mitosis_haploid_matrix(n, Model{u,k}, transition_t{});
                auto x = mitosis_matrix(n, Model{u,k}, transition_t{},1);
                auto expected = make_test_range(m);
                auto test = make_test_range(x);
                CHECK_CLOSE_RANGES(test, expected, prec);
            }
            for(int mutype : {0,1,2}) {
                BOOST_TEST_CONTEXT("mutype=" << mutype) {
                    auto m = mitosis_haploid_matrix(n, Model{u,k}, mutype);
                    auto x = mitosis_matrix(n, Model{u,k}, mutype,1);
                    auto expected = make_test_range(m);
                    auto test = make_test_range(x);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
            }
            BOOST_TEST_CONTEXT("mutype=mean") {
                auto m = mitosis_haploid_matrix(n, Model{u,k}, mean_t{});
                auto x = mitosis_matrix(n, Model{u,k}, mean_t{}, 1);
                auto expected = make_test_range(m);
                auto test = make_test_range(x);
                CHECK_CLOSE_RANGES(test, expected, prec);
            }
        }
        // Ploidy 2
        BOOST_TEST_CONTEXT("ploidy=2") {
            BOOST_TEST_CONTEXT("mutype=all") {
                auto m = mitosis_diploid_matrix(n, Model{u,k}, transition_t{});
                auto x = mitosis_matrix(n, Model{u,k}, transition_t{},2);
                auto expected = make_test_range(m);
                auto test = make_test_range(x);
                CHECK_CLOSE_RANGES(test, expected, prec);
            }
            for(int mutype : {0,1,2,3,4}) {
                BOOST_TEST_CONTEXT("mutype=" << mutype) {
                    auto m = mitosis_diploid_matrix(n, Model{u,k}, mutype);
                    auto x = mitosis_matrix(n, Model{u,k}, mutype, 2);
                    auto expected = make_test_range(m);
                    auto test = make_test_range(x);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
            }
            BOOST_TEST_CONTEXT("mutype=mean") {
                auto m = mitosis_diploid_matrix(n, Model{u,k}, mean_t{});
                auto x = mitosis_matrix(n, Model{u,k}, mean_t{}, 2);
                auto expected = make_test_range(m);
                auto test = make_test_range(x);
                CHECK_CLOSE_RANGES(test, expected, prec);
            }           
        }
    }};

    run_mutation_tests(test);
}


BOOST_AUTO_TEST_CASE(test_gamete_matrix) {
    auto test = [](int n, double u, double k, double prec) -> void {
        BOOST_TEST_CONTEXT("n=" << n << ", u=" << u << ", k=" << k)
    {
        auto f = Model{u,k}.TransitionMatrix(n);
        int n2 = (n*(n+1)/2);

        // Ploidy 1
        BOOST_TEST_CONTEXT("ploidy=1") {
            BOOST_TEST_CONTEXT("mutype=all") {
                auto m = mitosis_haploid_matrix(n, Model{u,k}, transition_t{});
                auto x = gamete_matrix(n, Model{u,k}, transition_t{}, 1);
                auto expected = make_test_range(m);
                auto test = make_test_range(x);
                CHECK_CLOSE_RANGES(test, expected, prec);
            }
            BOOST_TEST_CONTEXT("mutype=mean") {
                auto m = mitosis_haploid_matrix(n, Model{u,k}, mean_t{});
                auto x = gamete_matrix(n, Model{u,k}, mean_t{}, 1);
                auto expected = make_test_range(m);
                auto test = make_test_range(x);
                CHECK_CLOSE_RANGES(test, expected, prec);
            }
            for(int mutype : {0,1,2}) {
                BOOST_TEST_CONTEXT("mutype=" << mutype) {
                    auto m = mitosis_haploid_matrix(n, Model{u,k}, mutype);
                    auto x = gamete_matrix(n, Model{u,k}, mutype, 1);
                    auto expected = make_test_range(m);
                    auto test = make_test_range(x);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
            }
        }
        // Ploidy 2
        BOOST_TEST_CONTEXT("ploidy=2") {
            BOOST_TEST_CONTEXT("mutype=all") {
                auto m = meiosis_haploid_matrix(n, Model{u,k}, transition_t{});
                auto x = gamete_matrix(n, Model{u,k}, transition_t{}, 2);
                auto expected = make_test_range(m);
                auto test = make_test_range(x);
                CHECK_CLOSE_RANGES(test, expected, prec);
            }
            BOOST_TEST_CONTEXT("mutype=mean") {
                auto m = meiosis_haploid_matrix(n, Model{u,k}, mean_t{});
                auto x = gamete_matrix(n, Model{u,k}, mean_t{}, 2);
                auto expected = make_test_range(m);
                auto test = make_test_range(x);
                CHECK_CLOSE_RANGES(test, expected, prec);
            }
            for(int mutype : {0,1,2}) {
                BOOST_TEST_CONTEXT("mutype=" << mutype) {
                    auto m = meiosis_haploid_matrix(n, Model{u,k}, mutype);
                    auto x = gamete_matrix(n, Model{u,k}, mutype, 2);
                    auto expected = make_test_range(m);
                    auto test = make_test_range(x);
                    CHECK_CLOSE_RANGES(test, expected, prec);
                }
            }
        }
    }};

    run_mutation_tests(test);
}



BOOST_AUTO_TEST_CASE(test_number_of_parent_genotypes) {
    BOOST_CHECK_EQUAL(number_of_parent_genotypes(4, 1), 4);
    BOOST_CHECK_EQUAL(number_of_parent_genotypes(4, 2), 10);

    BOOST_CHECK_EQUAL(number_of_parent_genotype_pairs(4,1,1), 16);
    BOOST_CHECK_EQUAL(number_of_parent_genotype_pairs(4,1,2), 40);
    BOOST_CHECK_EQUAL(number_of_parent_genotype_pairs(4,2,1), 40);
    BOOST_CHECK_EQUAL(number_of_parent_genotype_pairs(4,2,2), 100);    
}

BOOST_AUTO_TEST_CASE(test_meiosis_matrix) {
    auto test_ploidy = [](int n, const int dad_ploidy, Model dad_m,
    const int mom_ploidy, Model mom_m, double prec) -> void {
        BOOST_TEST_CONTEXT("dad_ploidy=" << dad_ploidy << ", mom_ploidy=" << mom_ploidy)
    {
        auto dad_g = gamete_matrix(n, dad_m, transition_t{}, dad_ploidy);
        auto mom_g = gamete_matrix(n, mom_m, transition_t{}, mom_ploidy);

        // All Mutations
        std::vector<double> m;
        for(int x=0;x<n;++x) {
             for(int y=0;y<=x;++y) {
                for(int a=0;a<number_of_parent_genotypes(n,dad_ploidy);++a) {
                    for(int b=0;b<number_of_parent_genotypes(n,mom_ploidy);++b) {
                        // {a,b} => {x,y}
                        double d = 0.0;
                        d += dad_g(a,x)*mom_g(b,y);
                        if(x != y) {
                            d += dad_g(a,y)*mom_g(b,x);
                        }
                        m.push_back(d);
                    }
                }
            }
        }
        auto x_all = meiosis_matrix(n, dad_m, mom_m, transition_t{}, dad_ploidy, mom_ploidy);
        auto expected_all = m;
        auto test_all = make_test_range(x_all);
        CHECK_CLOSE_RANGES(test_all, expected_all, prec);
        
        // Mean Mutations
        auto x_mean = meiosis_matrix(n, dad_m, mom_m, mean_t{}, dad_ploidy, mom_ploidy);
        auto dad_a = gamete_matrix(n, dad_m, mean_t{}, dad_ploidy);
        auto mom_a = gamete_matrix(n, mom_m, mean_t{}, mom_ploidy);
        for(int x=0,z=0;x<n;++x) {
             for(int y=0;y<=x;++y) {
                for(int a=0;a<number_of_parent_genotypes(n,dad_ploidy);++a) {
                    for(int b=0;b<number_of_parent_genotypes(n,mom_ploidy);++b) {
                        // {a,b} => {x,y}
                        double d = 0.0;
                        d += dad_a(a,x)*mom_g(b,y) + dad_g(a,x)*mom_a(b,y);
                        if(x != y) {
                            d += dad_a(a,y)*mom_g(b,x) + dad_g(a,y)*mom_a(b,x);
                        }
                        m[z++] = d;
                    }
                }
            }
        }
        auto expected_mean = m;
        auto test_mean = make_test_range(x_mean);
        CHECK_CLOSE_RANGES(test_mean, expected_mean, prec);

        // Check for number of mutations
        for(int mutype : {0,1,2,3}) {
        BOOST_TEST_CONTEXT("mutype=" << mutype) {
            auto x_mutype = meiosis_matrix(n, dad_m, mom_m, mutype, dad_ploidy, mom_ploidy);
            m.assign(m.size(),0.0);
            for(int u=0;u<=mutype;++u) {
                auto dad_u = gamete_matrix(n, dad_m, mutype, dad_ploidy);
                auto mom_u = gamete_matrix(n, mom_m, mutype-u, mom_ploidy);
                for(int x=0,z=0;x<n;++x) {
                     for(int y=0;y<=x;++y) {
                        for(int a=0;a<number_of_parent_genotypes(n,dad_ploidy);++a) {
                            for(int b=0;b<number_of_parent_genotypes(n,mom_ploidy);++b) {
                                // {a,b} => {x,y}
                                double d = 0.0;
                                d += dad_u(a,x)*mom_u(b,y);
                                if(x != y) {
                                    d += dad_u(a,y)*mom_u(b,x);
                                }
                                m[z++] += d;
                            }
                        }
                    }
                }
            }
            auto expected_mutype = m;
            auto test_mutype = make_test_range(x_mutype);
            CHECK_CLOSE_RANGES(test_mutype, expected_mutype, prec);            
        }}
    }};

    auto test = [&](int n, double dad_u, double mom_u, double k, double prec) -> void {
        BOOST_TEST_CONTEXT("n=" << n << ", dad_u=" << dad_u
                                <<  ", mom_u=" << mom_u << ", k=" << k)
    {
        auto dad = Model{dad_u, k};
        auto mom = Model{mom_u, k};

        test_ploidy(n, 2, dad, 2, mom, prec);
        test_ploidy(n, 1, dad, 2, mom, prec);
        test_ploidy(n, 1, dad, 1, mom, prec);
        test_ploidy(n, 2, dad, 1, mom, prec);
    }};

    double prec = 2*DBL_EPSILON;
    test(4, 0.0, 0.0, 4, prec);

    test(1, 1e-8, 1e-7, 3, prec);
    test(2, 1e-7, 1e-8, 4, prec);
    test(3,  0.0, 1e-3, 5, prec);
    test(4, 1e-6, 1e-9, 6, prec);

}

BOOST_AUTO_TEST_CASE(test_population_prior_check) {
    BOOST_CHECK_EQUAL(population_prior_check(0.001, 0.0, 0.0, 0.0, 5.0), true);
    BOOST_CHECK_EQUAL(population_prior_check(0.0, 1.0, 1.0, 1.0, 5.0), true);
    BOOST_CHECK_EQUAL(population_prior_check(0.1, -20.0, -20.0, -10.0, 5.0), true);
    BOOST_CHECK_EQUAL(population_prior_check(100, 0.0, 0.0, 0.0, 5.0), true);

    BOOST_CHECK_EQUAL(population_prior_check(-0.1, 0.0, 0.0, 0.0, 5.0), false);
    BOOST_CHECK_EQUAL(population_prior_check(0.1, 2.0, 0.0, 0.0, 5.0), false);
    BOOST_CHECK_EQUAL(population_prior_check(0.1, 0.0, 2.0, 0.0, 5.0), false);
    BOOST_CHECK_EQUAL(population_prior_check(0.1, 0.0, 0.0, 2.0, 5.0), false);
    BOOST_CHECK_EQUAL(population_prior_check(0.1, -100.0, 0.0, 0.0, 5.0), false);
    BOOST_CHECK_EQUAL(population_prior_check(0.1, 0.0, -100.0, 0.0, 5.0), false);
    BOOST_CHECK_EQUAL(population_prior_check(0.1, 0.0, 0.0, -100.0, 5.0), false);
    BOOST_CHECK_EQUAL(population_prior_check(0.1, 0.0, 0.0, 0.0, 1.1), false);
}


template<typename F>
void run_population_tests(F test, double prec = 2*DBL_EPSILON) {
    auto test_prior = [&](d4 p) {
        test(0.001,  {0.25, 0.25, 0.25, 0.25}, p, prec);
        test(0.1, {0.25, 0.25, 0.25, 0.25}, p, prec);
        test(1e-6, {0.3,0.2,0.2,0.3}, p, prec);
        test(0.001, {0.1, 0.2, 0.3, 0.4}, p, prec); 
        test(0.001, {0.01, 0.1, 0.19, 0.7}, p, prec);
        test(1, {0.01, 0.1, 0.19, 0.7}, p, prec);
        test(100, {0.01, 0.1, 0.19, 0.7}, p, prec);        
        test(0.01,  {1, 1, 1, 1}, p, prec);
        test(0.0,  {1, 1, 1, 1}, p, prec);
        test(1e-9,  {4, 3, 2, 1}, p, prec);
    };
    test_prior({0.0,0,0,0});
    for(double w : {1.0,1e5,1e-5}) {
        test_prior({w,0,0,0});
        test_prior({0,w,0,0});
        test_prior({0,0,w,0});
        test_prior({0,0,0,w});
    }
}

BOOST_AUTO_TEST_CASE(test_population_prior_diploid) {
    constexpr double prec = 2*DBL_EPSILON;

    auto test_diploid = [prec](int num_obs_alleles, double theta, double hom_bias, double het_bias,
        double k_alleles) {
    BOOST_TEST_CONTEXT("num_obs_alleles=" << num_obs_alleles
        << ", theta=" << theta << ", hom_bias=" << hom_bias
        << ", het_bias=" << het_bias << ", k_alleles=" << k_alleles ) 
    {
        BOOST_REQUIRE(population_prior_check(theta, hom_bias, het_bias, 0.0, k_alleles));

        auto x_prior = population_prior_diploid(num_obs_alleles, theta, hom_bias,
            het_bias, k_alleles);
        if(num_obs_alleles == k_alleles) {
            BOOST_CHECK_CLOSE_FRACTION(x_prior.sum(), 1.0, prec);
        }
        std::vector<double> expected_prior;
        const int sz = num_obs_alleles;
        for(int a=0;a<sz;++a) {
            for(int b=0;b<=a;++b) {
                double k = k_alleles;
                double e = theta/(k_alleles-1.0);

                double r = 0.0;
                if(a == 0 && b == 0) {
                    // Reference Homozygote
                    r = (1.0+e)/(1.0+k*e)*(2.0+e+(k-1.0)*e*hom_bias)/(2.0+k*e);
                } else if(a == b) {
                    r = (1.0+e)/(1.0+k*e)*(e-e*hom_bias)/(2.0+k*e);
                } else if(b == 0 || a == 0) {
                    // Reference Heterozygote
                    r = e/(1.0+k*e)*(2.0+2.0*e+(k-2.0)*e*het_bias)/(2.0+k*e);
                } else {
                    // Alt heterozygote
                    r = e/(1.0+k*e)*(2.0*e-2.0*e*het_bias)/(2.0+k*e);
                }
                expected_prior.push_back(r);
            }
        }
        auto test_prior = make_test_range(x_prior);
        CHECK_CLOSE_RANGES(test_prior, expected_prior, prec);        
    }};

    auto test = [&](double theta, double hom_bias, double het_bias, double k_alleles) {
        for(int i=1;i<8;++i) {
            test_diploid(i, theta, hom_bias, het_bias, k_alleles);
        }
    }; 

    test(0.001, 0.0, 0.0, 5.0);
    test(0.01, 1.0, 0.0, 6.0);
    test(0.1, 0.0, 1.0, 4.0);
    test(0.001, 1.0, 1.0, 5.0);
    test(0.001, -1.0, -1.0, 5.5);
}

BOOST_AUTO_TEST_CASE(test_population_prior_haploid) {
    constexpr double prec = 2*DBL_EPSILON;

    auto test_haploid = [prec](int num_obs_alleles, double theta, double hap_bias,
     double k_alleles) {
    BOOST_TEST_CONTEXT("num_obs_alleles=" << num_obs_alleles
        << ", theta=" << theta << ", hap_bias=" << hap_bias
        << ", k_alleles=" << k_alleles ) 
    {
        BOOST_REQUIRE(population_prior_check(theta, 0.0, 0.0, hap_bias, k_alleles));

        auto x_prior = population_prior_haploid(num_obs_alleles, theta, hap_bias, k_alleles);
        if(num_obs_alleles == k_alleles) {
            BOOST_CHECK_CLOSE_FRACTION(x_prior.sum(), 1.0, prec);
        }

        std::vector<double> expected_prior;
        const int sz = num_obs_alleles;
        for(int a=0;a<sz;++a) {
            double k = k_alleles;
            double e = theta/(k_alleles-1.0);

            double r = 0.0;
            if(a == 0) {
                // reference haploid
                r = (1.0+e+e*(k-1)*hap_bias)/(1.0+k*e);
            } else {
                r = (e-e*hap_bias)/(1.0+k*e);
            }
            expected_prior.push_back(r);
        }
        auto test_prior = make_test_range(x_prior);
        CHECK_CLOSE_RANGES(test_prior, expected_prior, prec);        
    }};

    auto test = [&](double theta, double hap_bias, double k_alleles) {
        for(int i=1;i<8;++i) {
            test_haploid(i, theta, hap_bias, k_alleles);
        }
    }; 

    test(0.001, 0.0, 5.0);
    test(0.01, 1.0, 6.0);
    test(0.1, 0.0, 4.0);
    test(0.001, 1.0, 5.0);
    test(0.001, -1.0, 5.5);
}
