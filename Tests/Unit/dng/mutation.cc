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

    test(1, 1e-8, 4, prec);
    test(2, 1e-8, 4, prec);
    test(3, 1e-8, 4, prec);
    test(4, 1e-8, 4, prec);
    
    test(4, 1e-8, 4, prec);
    test(4, 1e-9, 5, prec);
    test(4, 1e-6, 6, prec);
    test(4, 1e-3, 7, prec);
}

BOOST_AUTO_TEST_CASE(test_model) {
    auto test = [](int n, double u, double k, double prec) -> void {
        BOOST_TEST_CONTEXT("n=" << n << ", u=" << u << ", k=" << k)
    {
        BOOST_REQUIRE(n <= k);

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

        BOOST_TEST_CONTEXT("model::TransitionMatrix") {
            std::vector<double> expected_matrix;
            for(int i=0;i<n;++i) {
                for(int j=0;j<n;++j) {
                    expected_matrix.push_back(P(j,i));
                }
            }

            auto model = Model{u,k};
            auto test_matrix_ = model.TransitionMatrix(n);
            auto test_matrix = make_test_range(test_matrix_);
            CHECK_CLOSE_RANGES(test_matrix, expected_matrix, prec);
        }

        // Test Event Matrix
        m = identity_matrix<double>(n+1);
        mat J = Q+m;
        mat S;

        for(int x=0; x<=10;++x) {
        BOOST_TEST_CONTEXT("model::EventTransitionMatrix w/ x=" << x) {
            if(x == 0) {
                P = m*exp(-u);
                S = P*x;
                f = 1.0;
                t = 1.0;
            } else {
                mat a = prec_prod(m,J);
                m = a;
                t *= u;
                f *= x;                
                P = m*(exp(-u)*t/f);
                S += P*x;
            }
            std::vector<double> expected_matrix;
            for(int i=0;i<n;++i) {
                for(int j=0;j<n;++j) {
                    expected_matrix.push_back(P(j,i));
                }
            }

            auto model = Model{u,k};
            auto test_matrix_ = model.EventTransitionMatrix(n,x);
            auto test_matrix = make_test_range(test_matrix_);
            CHECK_CLOSE_RANGES(test_matrix, expected_matrix, 256*prec);
        }}

        BOOST_TEST_CONTEXT("model::MeanTransitionMatrix") {
            std::vector<double> expected_matrix;
            for(int i=0;i<n;++i) {
                for(int j=0;j<n;++j) {
                    expected_matrix.push_back(S(j,i));
                }
            }

            auto model = Model{u,k};
            auto test_matrix_ = model.MeanTransitionMatrix(n);
            auto test_matrix = make_test_range(test_matrix_);
            CHECK_CLOSE_RANGES(test_matrix, expected_matrix, 4096*prec);
        }
    }};

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
        auto m_0 = Model{u,k}.EventTransitionMatrix(n,0);
        auto expected_0 = make_test_range(m_0);
        auto test_0 = make_test_range(x_0);
        CHECK_CLOSE_RANGES(test_0, expected_0, prec);

        // Test 1 Mutation
        auto x_1 = mitosis_haploid_matrix(n, Model{u,k}, 1);
        auto m_1 = Model{u,k}.EventTransitionMatrix(n,1);
        auto expected_1 = make_test_range(m_1);
        auto test_1 = make_test_range(x_1);
        CHECK_CLOSE_RANGES(test_1, expected_1, prec);

        // Test 2 Mutations
        auto x_2 = mitosis_haploid_matrix(n, Model{u,k}, 2);
        auto m_2 = Model{u,k}.EventTransitionMatrix(n,2);
        auto expected_2 = make_test_range(m_2);
        auto test_2 = make_test_range(x_2);
        CHECK_CLOSE_RANGES(test_2, expected_2, prec);

        // Test Mean Mutations
        auto x_mean = mitosis_haploid_matrix(n, Model{u,k}, mean_t{});
        auto m_mean = Model{u,k}.MeanTransitionMatrix(n);
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
        auto ff = Model{u,k}.TransitionMatrix(n);
        auto fa = Model{u,k}.MeanTransitionMatrix(n);
        auto f0 = Model{u,k}.EventTransitionMatrix(n, 0);
        auto f1 = Model{u,k}.EventTransitionMatrix(n, 1);
        auto f2 = Model{u,k}.EventTransitionMatrix(n, 2);

        // compute n*n x n*n mutation matrix
        boost::multi_array<double, 4> mf(boost::extents[n][n][n][n]);
        boost::multi_array<double, 4> ma(boost::extents[n][n][n][n]);
        boost::multi_array<double, 4> m0(boost::extents[n][n][n][n]);
        boost::multi_array<double, 4> m1(boost::extents[n][n][n][n]);
        boost::multi_array<double, 4> m2(boost::extents[n][n][n][n]);

        for(int a=0;a<n;++a) {
            for(int b=0;b<n;++b) {
                for(int x=0;x<n;++x) {
                    for(int y=0;y<n;++y) {
                        // a -> x and b -> y
                        mf[a][b][x][y] = ff(a,x)*ff(b,y);
                        ma[a][b][x][y] = ff(a,x)*fa(b,y)+fa(a,x)*ff(b,y);
                        m0[a][b][x][y] = f0(a,x)*f0(b,y);
                        m1[a][b][x][y] = f0(a,x)*f1(b,y)+f1(a,x)*f0(b,y);
                        m2[a][b][x][y] = f0(a,x)*f2(b,y)+f1(a,x)*f1(b,y)+f2(a,x)*f0(b,y);
                    }
                }
            }
        }
        auto fold = [&](const boost::multi_array<double, 4>& m) -> std::vector<double> {
            int n2 = (n*(n+1)/2);
            std::vector<double> r( n2*n2 );
            for(int x=0,z=0;x<n;++x) {
                for(int y=0;y<=x;++y) {
                    for(int a=0,c=0;a<n;++a) {
                        for(int b=0;b<=a;++b) {
                            r[z] = m[a][b][x][y];
                            if( x != y) {
                                r[z] += m[a][b][y][x];
                            }
                            ++z;
                        }
                    }
                }
            }
            return r;
        };

        // Test MUTATIONS_ALL
        auto x_all = mitosis_diploid_matrix(n, Model{u,k}, transition_t{});
        auto expected_all = fold(mf);
        auto test_all = make_test_range(x_all);
        CHECK_CLOSE_RANGES(test_all, expected_all, prec);

        // Test 0 Mutations
        auto x_0 = mitosis_diploid_matrix(n, Model{u,k}, 0);
        auto expected_0 = fold(m0);
        auto test_0 = make_test_range(x_0);
        CHECK_CLOSE_RANGES(test_0, expected_0, prec);

        // Test 1 Mutations
        auto x_1 = mitosis_diploid_matrix(n, Model{u,k}, 1);
        auto expected_1 = fold(m1);
        auto test_1 = make_test_range(x_1);
        CHECK_CLOSE_RANGES(test_1, expected_1, prec);

        // Test 2 Mutations
        auto x_2 = mitosis_diploid_matrix(n, Model{u,k}, 2);
        auto expected_2 = fold(m2);
        auto test_2 = make_test_range(x_2);
        CHECK_CLOSE_RANGES(test_2, expected_2, prec);

        // Test Mean Mutations
        auto x_mean = mitosis_diploid_matrix(n, Model{u,k}, mean_t{});
        auto expected_mean = fold(ma);
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
        auto ff = Model{u,k}.TransitionMatrix(n);
        auto fa = Model{u,k}.MeanTransitionMatrix(n);
        auto f0 = Model{u,k}.EventTransitionMatrix(n, 0);
        auto f1 = Model{u,k}.EventTransitionMatrix(n, 1);
        auto f2 = Model{u,k}.EventTransitionMatrix(n, 2);

        auto fold = [&](const decltype(ff)& m) -> std::vector<double> {
            int n2 = (n*(n+1)/2);
            std::vector<double> r(n*n2);

            for(int x=0,z=0;x<n;++x) {
                for(int a=0,c=0;a<n;++a) {
                    for(int b=0;b<=a;++b) {
                        r[z] = 0.5*(m(a,x)+m(b,x));
                        ++z;
                    }
                }
            }
            return r;
        };

        // Test MUTATIONS_ALL
        auto x_all = meiosis_haploid_matrix(n, Model{u,k}, transition_t{});
        auto expected_all = fold(ff);
        auto test_all = make_test_range(x_all);
        CHECK_CLOSE_RANGES(test_all, expected_all, prec);

        // Test 0 Mutations
        auto x_0 = meiosis_haploid_matrix(n, Model{u,k}, 0);
        auto expected_0 = fold(f0);
        auto test_0 = make_test_range(x_0);
        CHECK_CLOSE_RANGES(test_0, expected_0, prec);

        // Test 1 Mutation
        auto x_1 = meiosis_haploid_matrix(n, Model{u,k}, 1);
        auto expected_1 = fold(f1);
        auto test_1 = make_test_range(x_1);
        CHECK_CLOSE_RANGES(test_1, expected_1, prec);

        // Test 2 Mutations
        auto x_2 = meiosis_haploid_matrix(n, Model{u,k}, 2);
        auto expected_2 = fold(f2);
        auto test_2 = make_test_range(x_2);
        CHECK_CLOSE_RANGES(test_2, expected_2, prec);

        // Test Mean Mutations
        auto x_mean = meiosis_haploid_matrix(n, Model{u,k}, mean_t{});
        auto expected_mean = fold(fa);
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
        int dad_sz = number_of_parent_genotypes(n,dad_ploidy);
        int mom_sz = number_of_parent_genotypes(n,mom_ploidy);

        auto dad_ff = gamete_matrix(n, dad_m, transition_t{}, dad_ploidy);
        auto dad_fa = gamete_matrix(n, dad_m, mean_t{}, dad_ploidy);
        auto dad_f0 = gamete_matrix(n, dad_m, 0, dad_ploidy);
        auto dad_f1 = gamete_matrix(n, dad_m, 1, dad_ploidy);
        auto dad_f2 = gamete_matrix(n, dad_m, 2, dad_ploidy);

        auto mom_ff = gamete_matrix(n, mom_m, transition_t{}, mom_ploidy);
        auto mom_fa = gamete_matrix(n, mom_m, mean_t{}, mom_ploidy);
        auto mom_f0 = gamete_matrix(n, mom_m, 0, mom_ploidy);
        auto mom_f1 = gamete_matrix(n, mom_m, 1, mom_ploidy);
        auto mom_f2 = gamete_matrix(n, mom_m, 2, mom_ploidy);

        auto ext = boost::extents[dad_sz][mom_sz][n][n];
        boost::multi_array<double, 4> mf(ext);
        boost::multi_array<double, 4> ma(ext);
        boost::multi_array<double, 4> m0(ext);
        boost::multi_array<double, 4> m1(ext);
        boost::multi_array<double, 4> m2(ext);

        for(int a=0;a<dad_sz;++a) {
            for(int b=0;b<mom_sz;++b) {
                for(int x=0;x<n;++x) {
                    for(int y=0;y<n;++y) {
                        // a -> x and b -> y
                        mf[a][b][x][y] = dad_ff(a,x)*mom_ff(b,y);
                        ma[a][b][x][y] = dad_ff(a,x)*mom_fa(b,y)+dad_fa(a,x)*mom_ff(b,y);
                        m0[a][b][x][y] = dad_f0(a,x)*mom_f0(b,y);
                        m1[a][b][x][y] = dad_f0(a,x)*mom_f1(b,y)+dad_f1(a,x)*mom_f0(b,y);
                        m2[a][b][x][y] = dad_f0(a,x)*mom_f2(b,y)+dad_f1(a,x)*mom_f1(b,y)+dad_f2(a,x)*mom_f0(b,y);
                    }
                }
            }
        }

        auto fold = [&](const boost::multi_array<double, 4>& m) -> std::vector<double> {
            int n2 = (n*(n+1)/2);
            std::vector<double> r( dad_sz*mom_sz*n2 );
            for(int x=0,z=0;x<n;++x) {
                for(int y=0;y<=x;++y) {
                    for(int a=0;a<dad_sz;++a) {
                        for(int b=0;b<mom_sz;++b) {
                            r[z] = m[a][b][x][y];
                            if( x != y) {
                                r[z] += m[a][b][y][x];
                            }
                            ++z;
                        }
                    }
                }
            }
            return r;
        };
        auto x_all = meiosis_matrix(n, dad_m, mom_m, transition_t{}, dad_ploidy, mom_ploidy);
        auto expected_all = fold(mf);
        auto test_all = make_test_range(x_all);
        CHECK_CLOSE_RANGES(test_all, expected_all, prec);
        
        // Mean Mutations
        auto x_mean = meiosis_matrix(n, dad_m, mom_m, mean_t{}, dad_ploidy, mom_ploidy);
        auto expected_mean = fold(ma);
        auto test_mean = make_test_range(x_mean);
        CHECK_CLOSE_RANGES(test_mean, expected_mean, prec);

        auto x_0 = meiosis_matrix(n, dad_m, mom_m, 0, dad_ploidy, mom_ploidy);
        auto expected_0 = fold(m0);
        auto test_0 = make_test_range(x_0);
        CHECK_CLOSE_RANGES(test_0, expected_0, prec);

        auto x_1 = meiosis_matrix(n, dad_m, mom_m, 1, dad_ploidy, mom_ploidy);
        auto expected_1 = fold(m1);
        auto test_1 = make_test_range(x_1);
        CHECK_CLOSE_RANGES(test_1, expected_1, prec);

        auto x_2 = meiosis_matrix(n, dad_m, mom_m, 2, dad_ploidy, mom_ploidy);
        auto expected_2 = fold(m2);
        auto test_2 = make_test_range(x_2);
        CHECK_CLOSE_RANGES(test_2, expected_2, prec);
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
