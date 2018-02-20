/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2016,2017 Reed A. Cartwright
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

#define BOOST_TEST_MODULE dng::peel

#include <dng/peeling.h>

#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

#include <boost/range/algorithm/generate.hpp>
#include <boost/range/algorithm/copy.hpp>

#include <dng/mutation.h>

#include "../testing.h"
#include "../xorshift64.h"

using namespace dng;
using namespace dng::peel;
using dng::detail::make_test_range;

const int NUM_TEST = 25;
int g_seed_counter = 0;

BOOST_AUTO_TEST_CASE(test_peel_workspace) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    workspace_t work;
    work.Resize(10);
    work.founder_nodes = {0,2};
    work.germline_nodes = {2,4};
    work.somatic_nodes = {4,5};
    work.library_nodes = {5,10};

    work.ploidies.assign(10,2);
    work.ploidies[1] = 1;
    work.ploidies[2] = 1;
    work.ploidies[4] = 1;

    // SetGermline
    GenotypeArray haploid_prior(4), diploid_prior(10);
    generate(make_test_range(haploid_prior), [&](){ return xrand.get_double52(); });
    generate(make_test_range(diploid_prior), [&](){ return xrand.get_double52(); });

    work.SetGermline(diploid_prior, haploid_prior);
    
    auto test_upper0 = make_test_range(work.upper[0]);
    auto expected_upper0 = make_test_range(diploid_prior);
    CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

    auto test_upper1 = make_test_range(work.upper[1]);
    auto expected_upper1 = make_test_range(haploid_prior);
    CHECK_EQUAL_RANGES(test_upper1, expected_upper1);

    // CleanupFast
    for(auto &&a : work.lower) {
        generate(make_test_range(a), [&](){ return xrand.get_double52(); });
    }
    work.dirty_lower = true;
    work.CleanupFast();
    BOOST_CHECK_EQUAL(work.dirty_lower, false);
    for(int i=0;i<work.somatic_nodes.second;++i) {
        BOOST_TEST_INFO("node=" << i);
        auto test_range = make_test_range(work.lower[i]);
        size_t sz = (i == 1 || i == 2) ? 4 : 10;
        std::vector<double> expected_range(sz,1); 
        CHECK_EQUAL_RANGES(test_range, expected_range);
    }
    for(int i=work.library_nodes.first;i<work.num_nodes;++i) {
        BOOST_TEST_INFO("node=" << i);
        auto test_range = make_test_range(work.lower[i]);
        std::vector<double> expected_range(10,1); 
        CHECK_NE_RANGES(test_range, expected_range);
    }

    // Cleanup
    for(auto &&a : work.lower) {
        generate(make_test_range(a), [&](){ return xrand.get_double52(); });
    }
    work.dirty_lower = true;
    work.Cleanup();
    BOOST_CHECK_EQUAL(work.dirty_lower, false);
    for(int i=0;i<work.num_nodes;++i) {
        BOOST_TEST_INFO("node=" << i);
        auto test_range = make_test_range(work.lower[i]);
        size_t sz = (i == 1 || i == 2) ? 4 : 10;
        std::vector<double> expected_range(sz,1); 
        CHECK_EQUAL_RANGES(test_range, expected_range);
    }
}

BOOST_AUTO_TEST_CASE(test_peel_up_fast) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto m = Mk::matrix(4, 1e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1),
                                expected_lower0(sz0);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(2);
            if(ploidy0 == ploidy1) {
                mats[1] = mitosis_matrix(ploidy0, m);
            } else {
                mats[1] = gamete_matrix(ploidy0, m);
            }

            // Calculated expected value
            for(int i=0;i<sz0;++i) {
                expected_lower0[i] = 0.0;
                for(int j=0;j<sz1;++j) {
                    expected_lower0[i] += mats[1](i,j) * lower1[j];
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(2);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());

            // do the peeling
            up_fast(workspace, {0,1}, mats);

            auto expected_lower1 = lower1;
            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_EQUAL_RANGES(test_lower1, expected_lower1);

            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_CLOSE_RANGES(test_lower0, expected_lower0, prec);        
        }}        
    };

    test(2,2);
    test(2,1);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_up) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto m = Mk::matrix(4, 1e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1),
                                expected_lower0(sz0);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(2);
            if(ploidy0 == ploidy1) {
                mats[1] = mitosis_matrix(ploidy0, m);
            } else {
                mats[1] = gamete_matrix(ploidy0, m);
            }

            // Calculated expected value
            for(int i=0;i<sz0;++i) {
                expected_lower0[i] = 0.0;
                for(int j=0;j<sz1;++j) {
                    expected_lower0[i] += mats[1](i,j) * lower1[j];
                }
                expected_lower0[i] *= lower0[i];
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(2);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());

            // do the peeling
            up(workspace, {0,1}, mats);

            auto expected_lower1 = lower1;
            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_EQUAL_RANGES(test_lower1, expected_lower1);

            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_CLOSE_RANGES(test_lower0, expected_lower0, prec);
        }}
    };

    test(2,2);
    test(2,1);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_down_fast) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto m = Mk::matrix(4, 1e-6, log(4));

            // setup data
            std::vector<double> upper0(sz0), upper1(sz1), expected_upper1(sz1);
            generate(upper0, [&](){ return xrand.get_double52(); });
            generate(upper1, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(2);
            if(ploidy0 == ploidy1) {
                mats[1] = mitosis_matrix(ploidy0, m);
            } else {
                mats[1] = gamete_matrix(ploidy0, m);
            }

            // Calculated expected value
            for(int i=0;i<sz1;++i) {
                expected_upper1[i] = 0.0;
                for(int j=0;j<sz0;++j) {
                    expected_upper1[i] += mats[1](j,i) * upper0[j];
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.upper.resize(2);
            workspace.upper[0].resize(sz0);
            workspace.upper[1].resize(sz1);
            copy(upper0, workspace.upper[0].data());
            copy(upper1, workspace.upper[1].data());

            // do the peeling
            down_fast(workspace, {0,1}, mats);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_CLOSE_RANGES(test_upper1, expected_upper1, prec);        
        }}
    };

    test(2,2);
    test(2,1);
    test(1,1);    
}

BOOST_AUTO_TEST_CASE(test_peel_down) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto m = Mk::matrix(4, 1e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), upper0(sz0), upper1(sz1), expected_upper1(sz1);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            generate(upper1, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(2);
            if(ploidy0 == ploidy1) {
                mats[1] = mitosis_matrix(ploidy0, m);
            } else {
                mats[1] = gamete_matrix(ploidy0, m);
            }

            // Calculated expected value
            for(int i=0;i<sz1;++i) {
                expected_upper1[i] = 0.0;
                for(int j=0;j<sz0;++j) {
                    expected_upper1[i] += mats[1](j,i) * upper0[j]*lower0[j];
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(1);
            workspace.lower[0].resize(sz0);
            copy(lower0, workspace.lower[0].data());
            workspace.upper.resize(2);
            workspace.upper[0].resize(sz0);
            workspace.upper[1].resize(sz1);
            copy(upper0, workspace.upper[0].data());
            copy(upper1, workspace.upper[1].data());

            // do the peeling
            down(workspace, {0,1}, mats);

            auto expected_lower0 = lower0;
            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_EQUAL_RANGES(test_lower0, expected_lower0);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_CLOSE_RANGES(test_upper1, expected_upper1, prec);        
        }}
    };

    test(2,2);
    test(2,1);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_tofather_fast) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        const int sz2 = 10;
        const int sz3 = 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto da = Mk::matrix(4, 1e-6, log(4));
            auto ma = Mk::matrix(4, 0.7e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), lower2(sz2), lower3(sz3), 
                    upper1(sz1), expected_lower0(sz0);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(lower2, [&](){ return xrand.get_double52(); });
            generate(lower3, [&](){ return xrand.get_double52(); });
            generate(upper1, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(4);
            mats[2] = meiosis_matrix(ploidy0,da,ploidy1,ma);
            mats[3] = meiosis_matrix(ploidy0,da,ploidy1,ma);

            // Calculated expected value
            for(int i=0;i<sz0;++i) {
                expected_lower0[i] = 0.0;
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    double d2 = 0.0;
                    for(int k=0;k<sz2;++k) {
                        d2 += mats[2](ij,k) * lower2[k];
                    }
                    double d3 = 0.0;
                    for(int k=0;k<sz2;++k) {
                        d3 += mats[3](ij,k) * lower3[k];
                    }
                    expected_lower0[i] += d2*d3*(lower1[j]*upper1[j]);
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(4);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.lower[2].resize(sz2);
            workspace.lower[3].resize(sz3);
            workspace.upper.resize(4);
            workspace.upper[1].resize(sz1);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(lower2, workspace.lower[2].data());
            copy(lower3, workspace.lower[3].data());
            copy(upper1, workspace.upper[1].data());

            // do the peeling
            to_father_fast(workspace, {0,1,2,3}, mats);

            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_CLOSE_RANGES(test_lower0, expected_lower0, prec);        

            auto expected_lower1 = lower1;
            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_EQUAL_RANGES(test_lower1, expected_lower1);

            auto expected_upper1 = upper1;
            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_EQUAL_RANGES(test_upper1, expected_upper1);

            auto expected_lower2 = lower2;
            auto test_lower2 = make_test_range(workspace.lower[2]);
            CHECK_EQUAL_RANGES(test_lower2, expected_lower2);

            auto expected_lower3 = lower3;
            auto test_lower3 = make_test_range(workspace.lower[3]);
            CHECK_EQUAL_RANGES(test_lower3, expected_lower3);
        }}
    };

    test(2,2);
    test(2,1);
    test(1,2);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_tofather) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        const int sz2 = 10;
        const int sz3 = 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto da = Mk::matrix(4, 1e-6, log(4));
            auto ma = Mk::matrix(4, 0.7e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), lower2(sz2), lower3(sz3), 
                    upper1(sz1), expected_lower0(sz0);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(lower2, [&](){ return xrand.get_double52(); });
            generate(lower3, [&](){ return xrand.get_double52(); });
            generate(upper1, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(4);
            mats[2] = meiosis_matrix(ploidy0,da,ploidy1,ma);
            mats[3] = meiosis_matrix(ploidy0,da,ploidy1,ma);

            // Calculated expected value
            for(int i=0;i<sz0;++i) {
                expected_lower0[i] = 0.0;
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    double d2 = 0.0;
                    for(int k=0;k<sz2;++k) {
                        d2 += mats[2](ij,k) * lower2[k];
                    }
                    double d3 = 0.0;
                    for(int k=0;k<sz2;++k) {
                        d3 += mats[3](ij,k) * lower3[k];
                    }
                    expected_lower0[i] += d2*d3*(lower1[j]*upper1[j]);
                }
                expected_lower0[i] *= lower0[i];
            }


            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(4);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.lower[2].resize(sz2);
            workspace.lower[3].resize(sz3);
            workspace.upper.resize(4);
            workspace.upper[1].resize(sz1);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(lower2, workspace.lower[2].data());
            copy(lower3, workspace.lower[3].data());
            copy(upper1, workspace.upper[1].data());

            // do the peeling
            to_father(workspace, {0,1,2,3}, mats);

            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_CLOSE_RANGES(test_lower0, expected_lower0, prec);        

            auto expected_lower1 = lower1;
            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_EQUAL_RANGES(test_lower1, expected_lower1);

            auto expected_upper1 = upper1;
            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_EQUAL_RANGES(test_upper1, expected_upper1);

            auto expected_lower2 = lower2;
            auto test_lower2 = make_test_range(workspace.lower[2]);
            CHECK_EQUAL_RANGES(test_lower2, expected_lower2);

            auto expected_lower3 = lower3;
            auto test_lower3 = make_test_range(workspace.lower[3]);
            CHECK_EQUAL_RANGES(test_lower3, expected_lower3);
        }}
    };

    test(2,2);
    test(2,1);
    test(1,2);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_tomother_fast) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        const int sz2 = 10;
        const int sz3 = 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto da = Mk::matrix(4, 1e-6, log(4));
            auto ma = Mk::matrix(4, 0.7e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), lower2(sz2), lower3(sz3), 
                    upper0(sz0), expected_lower1(sz1);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(lower2, [&](){ return xrand.get_double52(); });
            generate(lower3, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(4);
            mats[2] = meiosis_matrix(ploidy0,da,ploidy1,ma);
            mats[3] = meiosis_matrix(ploidy0,da,ploidy1,ma);

            // Calculated expected value
            for(int i=0;i<sz1;++i) {
                expected_lower1[i] = 0.0;
                for(int j=0;j<sz0;++j) {
                    int ij = j*sz1+i;
                    double d2 = 0.0;
                    for(int k=0;k<sz2;++k) {
                        d2 += mats[2](ij,k) * lower2[k];
                    }
                    double d3 = 0.0;
                    for(int k=0;k<sz2;++k) {
                        d3 += mats[3](ij,k) * lower3[k];
                    }
                    expected_lower1[i] += d2*d3*(lower0[j]*upper0[j]);
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(4);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.lower[2].resize(sz2);
            workspace.lower[3].resize(sz2);
            workspace.upper.resize(4);
            workspace.upper[0].resize(sz0);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(lower2, workspace.lower[2].data());
            copy(lower3, workspace.lower[3].data());
            copy(upper0, workspace.upper[0].data());

            // do the peeling
            to_mother_fast(workspace, {0,1,2,3}, mats);

            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_CLOSE_RANGES(test_lower1, expected_lower1, prec);        

            auto expected_lower0 = lower0;
            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_EQUAL_RANGES(test_lower0, expected_lower0);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto expected_lower2 = lower2;
            auto test_lower2 = make_test_range(workspace.lower[2]);
            CHECK_EQUAL_RANGES(test_lower2, expected_lower2);

            auto expected_lower3 = lower3;
            auto test_lower3 = make_test_range(workspace.lower[3]);
            CHECK_EQUAL_RANGES(test_lower3, expected_lower3);
        }}
    };

    test(2,2);
    test(2,1);
    test(1,2);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_tomother) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        const int sz2 = 10;
        const int sz3 = 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto da = Mk::matrix(4, 1e-6, log(4));
            auto ma = Mk::matrix(4, 0.7e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), lower2(sz2), lower3(sz3),
                    upper0(sz0), expected_lower1(sz1);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(lower2, [&](){ return xrand.get_double52(); });
            generate(lower3, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(4);
            mats[2] = meiosis_matrix(ploidy0,da,ploidy1,ma);
            mats[3] = meiosis_matrix(ploidy0,da,ploidy1,ma);

            // Calculated expected value
            for(int i=0;i<sz1;++i) {
                expected_lower1[i] = 0.0;
                for(int j=0;j<sz0;++j) {
                    int ij = j*sz1+i;
                    double d2 = 0.0;
                    for(int k=0;k<sz2;++k) {
                        d2 += mats[2](ij,k) * lower2[k];
                    }
                    double d3 = 0.0;
                    for(int k=0;k<sz2;++k) {
                        d3 += mats[3](ij,k) * lower3[k];
                    }
                    expected_lower1[i] += d2*d3*(lower0[j]*upper0[j]);
                }
                expected_lower1[i] *= lower1[i];
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(4);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.lower[2].resize(sz2);
            workspace.lower[3].resize(sz3);
            workspace.upper.resize(4);
            workspace.upper[0].resize(sz0);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(lower2, workspace.lower[2].data());
            copy(lower3, workspace.lower[3].data());
            copy(upper0, workspace.upper[0].data());

            // do the peeling
            to_mother(workspace, {0,1,2,3}, mats);

            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_CLOSE_RANGES(test_lower1, expected_lower1, prec);        

            auto expected_lower0 = lower0;
            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_EQUAL_RANGES(test_lower0, expected_lower0);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto expected_lower2 = lower2;
            auto test_lower2 = make_test_range(workspace.lower[2]);
            CHECK_EQUAL_RANGES(test_lower2, expected_lower2);

            auto expected_lower3 = lower3;
            auto test_lower3 = make_test_range(workspace.lower[3]);
            CHECK_EQUAL_RANGES(test_lower3, expected_lower3);
        }}
    };

    test(2,2);
    test(2,1);
    test(1,2);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_tochild_fast) {
    const double prec = 16.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        const int sz2 = 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto da = Mk::matrix(4, 1e-6, log(4));
            auto ma = Mk::matrix(4, 0.7e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), lower2(sz2), 
                                upper0(sz0), upper1(sz1), upper2(sz2),
                                expected_upper2(sz2);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(lower2, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            generate(upper1, [&](){ return xrand.get_double52(); });
            generate(upper2, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(3);
            mats[2] = meiosis_matrix(ploidy0,da,ploidy1,ma);

            // Calculated expected value
            for(int k=0;k<sz2;++k) {
                expected_upper2[k] = 0.0;
                for(int i=0,ij=0;i<sz0;++i) {
                    for(int j=0;j<sz1;++j,++ij) {
                        expected_upper2[k] += mats[2](ij,k) * (lower0[i]*upper0[i]) * (lower1[j]*upper1[j]);
                    }
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(3);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.lower[2].resize(sz2);
            workspace.upper.resize(3);
            workspace.upper[0].resize(sz0);
            workspace.upper[1].resize(sz1);
            workspace.upper[2].resize(sz2);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(lower2, workspace.lower[2].data());
            copy(upper0, workspace.upper[0].data());
            copy(upper1, workspace.upper[1].data());

            // do the peeling
            to_child_fast(workspace, {0,1,2}, mats);

            auto test_upper2 = make_test_range(workspace.upper[2]);
            CHECK_CLOSE_RANGES(test_upper2, expected_upper2, prec);

            auto expected_lower0 = lower0;
            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_EQUAL_RANGES(test_lower0, expected_lower0);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto expected_lower1 = lower1;
            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_EQUAL_RANGES(test_lower1, expected_lower1);

            auto expected_upper1 = upper1;
            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_EQUAL_RANGES(test_upper1, expected_upper1);

            auto expected_lower2 = lower2;
            auto test_lower2 = make_test_range(workspace.lower[2]);
            CHECK_EQUAL_RANGES(test_lower2, expected_lower2);
        }}
    };

    test(2,2);
    test(2,1);
    test(1,2);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_tochild) {
    const double prec = 32.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        const int sz2 = 10;
        const int sz3 = 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto da = Mk::matrix(4, 1e-6, log(4));
            auto ma = Mk::matrix(4, 0.7e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), lower2(sz2), lower3(sz3), 
                                upper0(sz0), upper1(sz1), upper2(sz2), upper3(sz3),
                                expected_upper2(sz2);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(lower2, [&](){ return xrand.get_double52(); });
            generate(lower3, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            generate(upper1, [&](){ return xrand.get_double52(); });
            generate(upper2, [&](){ return xrand.get_double52(); });
            generate(upper3, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(4);
            mats[2] = meiosis_matrix(ploidy0,da,ploidy1,ma);
            mats[3] = meiosis_matrix(ploidy0,ma,ploidy1,da);

            // Calculated expected value
            for(int k=0;k<sz2;++k) {
                expected_upper2[k] = 0.0;
                for(int i=0,ij=0;i<sz0;++i) {
                    for(int j=0;j<sz1;++j,++ij) {
                        for(int h=0;h<sz3;++h) {
                            expected_upper2[k] += mats[2](ij,k)
                            * (lower0[i]*upper0[i]) * (lower1[j]*upper1[j])
                            * (mats[3](ij,h)*lower3[h]);
                        }
                    }
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(4);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.lower[2].resize(sz2);
            workspace.lower[3].resize(sz3);
            workspace.upper.resize(4);
            workspace.upper[0].resize(sz0);
            workspace.upper[1].resize(sz1);
            workspace.upper[2].resize(sz2);
            workspace.upper[3].resize(sz3);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(lower2, workspace.lower[2].data());
            copy(lower3, workspace.lower[3].data());
            copy(upper0, workspace.upper[0].data());
            copy(upper1, workspace.upper[1].data());

            // do the peeling
            to_child(workspace, {0,1,2,3}, mats);

            auto test_upper2 = make_test_range(workspace.upper[2]);
            CHECK_CLOSE_RANGES(test_upper2, expected_upper2, prec);

            auto expected_lower0 = lower0;
            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_EQUAL_RANGES(test_lower0, expected_lower0);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto expected_lower1 = lower1;
            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_EQUAL_RANGES(test_lower1, expected_lower1);

            auto expected_upper1 = upper1;
            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_EQUAL_RANGES(test_upper1, expected_upper1);

            auto expected_lower2 = lower2;
            auto test_lower2 = make_test_range(workspace.lower[2]);
            CHECK_EQUAL_RANGES(test_lower2, expected_lower2);

            auto expected_lower3 = lower3;
            auto test_lower3 = make_test_range(workspace.lower[3]);
            CHECK_EQUAL_RANGES(test_lower3, expected_lower3);
        }}
    };

    test(2,2);
    test(2,1);
    test(1,2);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_up_reverse) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto m =  Mk::matrix(4, 1e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), upper0(sz0),
                                expected_super1(sz0), expected_upper1(sz1),
                                expected_lower0(sz0);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(2);
            if(ploidy0 == ploidy1) {
                mats[1] = mitosis_matrix(ploidy0, m);
            } else {
                mats[1] = gamete_matrix(ploidy0, m);
            }

            // Calculated expected value
            for(int i=0;i<sz0;++i) {
                expected_super1[i] = upper0[i]*lower0[i];
            }
            for(int j=0;j<sz1;++j) {
                expected_upper1[j] = 0.0;
                for(int i=0;i<sz0;++i) {
                    expected_upper1[j] += mats[1](i,j)*expected_super1[i];
                }
            }           

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(2);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.upper.resize(2);
            workspace.upper[0].resize(sz0);
            workspace.upper[1].resize(sz1);
            workspace.super.resize(2);
            workspace.super[1].resize(sz0);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(upper0, workspace.upper[0].data());

            // do the peeling
            up(workspace, {0,1}, mats);
            copy_n(workspace.lower[0].data(), sz0, expected_lower0.begin());
            up_reverse(workspace, {0,1}, mats);

            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_EQUAL_RANGES(test_lower0, expected_lower0);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto expected_lower1 = lower1;
            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_EQUAL_RANGES(test_lower1, expected_lower1);

            auto test_super1 = make_test_range(workspace.super[1]);
            CHECK_CLOSE_RANGES(test_super1, expected_super1, prec);

            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_CLOSE_RANGES(test_upper1, expected_upper1, prec);
        }}
    };

    test(2,2);
    test(2,1);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_down_reverse) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", test_num=" << test_num) {
            auto m =  Mk::matrix(4, 1e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), upper0(sz0),
                                expected_super1(sz0), expected_lower0(sz0);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(2);
            if(ploidy0 == ploidy1) {
                mats[1] = mitosis_matrix(ploidy0, m);
            } else {
                mats[1] = gamete_matrix(ploidy0, m);
            }

            // Calculated expected value
            for(int i=0;i<sz0;++i) {
                expected_super1[i] = upper0[i]*lower0[i];   
            }
            for(int i=0;i<sz0;++i) {
                expected_lower0[i] = 0.0;
                for(int j=0;j<sz1;++j) {
                    expected_lower0[i] += mats[1](i,j)*lower1[j];
                }
                expected_lower0[i] *= lower0[i];
            }           

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(2);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.upper.resize(2);
            workspace.upper[0].resize(sz0);
            workspace.upper[1].resize(sz1);
            workspace.super.resize(2);
            workspace.super[1].resize(sz0);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(upper0, workspace.upper[0].data());

            // do the peeling
            down(workspace, {0,1}, mats);
            down_reverse(workspace, {0,1}, mats);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto expected_lower1 = lower1;
            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_EQUAL_RANGES(test_lower1, expected_lower1);

            auto test_super1 = make_test_range(workspace.super[1]);
            CHECK_CLOSE_RANGES(test_super1, expected_super1, prec);

            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_CLOSE_RANGES(test_lower0, expected_lower0, prec);
        }}
    };

    test(2,2);
    test(2,1);
    test(1,1);
}

BOOST_AUTO_TEST_CASE(test_peel_tofather_reverse) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1, int ploidy2, int ploidy3) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        const int sz2 = (ploidy2 == 1) ? 4 : 10;
        const int sz3 = (ploidy3 == 1) ? 4 : 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", ploidy2=" << ploidy2
                        << ", ploidy3=" << ploidy3
                        << ", test_num=" << test_num) {
            auto da = Mk::matrix(4, 1e-6, log(4));
            auto ma = Mk::matrix(4, 0.7e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), lower2(sz2), lower3(sz3), 
                                upper0(sz0), upper1(sz1),
                                expected_lower0(sz0), expected_lower1(sz1),
                                expected_super2(sz0*sz1), expected_upper2(sz2),
                                expected_super3(sz0*sz1), expected_upper3(sz3);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(lower2, [&](){ return xrand.get_double52(); });
            generate(lower3, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            generate(upper1, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(4);
            mats[2] = meiosis_matrix(ploidy0,da,ploidy1,ma);
            mats[3] = meiosis_matrix(ploidy0,ma,ploidy1,da);

            std::vector<double> p(sz0*sz1);

            // Calculate expected value
            for(int i=0;i<sz0;++i) {
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    p[ij] = 1.0;
                    double d = 0.0;
                    for(int k=0;k<sz2;++k) {           
                        d += mats[2](ij,k) * lower2[k];
                    }
                    p[ij] *= d;
                    d = 0.0;
                    for(int k=0;k<sz3;++k) {           
                        d += mats[3](ij,k) * lower3[k];
                    }
                    p[ij] *= d;
                }
            }

            for(int j=0;j<sz1;++j) {
                expected_lower1[j] = 0.0;
                for(int i=0;i<sz0;++i) {
                    int ij = i*sz1+j;
                    expected_lower1[j] += p[ij]*(upper0[i]*lower0[i]);
                }
                expected_lower1[j] *= lower1[j];
            }

            for(int i=0;i<sz0;++i) {
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    double d = 0.0;
                    for(int k=0;k<sz3;++k) {
                        d += mats[3](ij,k) * lower3[k];
                    }
                    expected_super2[ij] = d
                        *(upper0[i]*lower0[i])
                        *(lower1[j]*upper1[j]);
                }
            }

            for(int i=0;i<sz0;++i) {
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    double d = 0.0;
                    for(int k=0;k<sz2;++k) {           
                        d += mats[2](ij,k) * lower2[k];
                    }
                    expected_super3[ij] = d
                        *(upper0[i]*lower0[i])
                        *(lower1[j]*upper1[j]);
                }
            }

            for(int k=0;k<sz2;++k) {
                expected_upper2[k] = 0.0;
                for(int ij=0;ij<sz0*sz1;++ij) {
                    expected_upper2[k] += expected_super2[ij]*mats[2](ij,k);
                }
            }

            for(int k=0;k<sz3;++k) {
                expected_upper3[k] = 0.0;
                for(int ij=0;ij<sz0*sz1;++ij) {
                    expected_upper3[k] += expected_super3[ij]*mats[3](ij,k);
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(4);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.lower[2].resize(sz2);
            workspace.lower[3].resize(sz3);
            workspace.upper.resize(4);
            workspace.upper[0].resize(sz0);
            workspace.upper[1].resize(sz1);
            workspace.super.resize(4);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(lower2, workspace.lower[2].data());
            copy(lower3, workspace.lower[3].data());
            copy(upper0, workspace.upper[0].data());
            copy(upper1, workspace.upper[1].data());

            // do the peeling
            to_father(workspace, {0,1,2,3}, mats);
            copy_n(workspace.lower[0].data(), sz0, expected_lower0.begin());
            to_father_reverse(workspace, {0,1,2,3}, mats);

            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_CLOSE_RANGES(test_lower1, expected_lower1, prec);

            auto test_super2 = make_test_range(workspace.super[2]);
            CHECK_CLOSE_RANGES(test_super2, expected_super2, prec);

            auto test_super3 = make_test_range(workspace.super[3]);
            CHECK_CLOSE_RANGES(test_super3, expected_super3, prec);

            auto test_upper2 = make_test_range(workspace.upper[2]);
            CHECK_CLOSE_RANGES(test_upper2, expected_upper2, prec);

            auto test_upper3 = make_test_range(workspace.upper[3]);
            CHECK_CLOSE_RANGES(test_upper3, expected_upper3, prec);

            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_EQUAL_RANGES(test_lower0, expected_lower0);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto expected_upper1 = upper1;
            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_EQUAL_RANGES(test_upper1, expected_upper1);

            auto expected_lower2 = lower2;
            auto test_lower2 = make_test_range(workspace.lower[2]);
            CHECK_EQUAL_RANGES(test_lower2, expected_lower2);

            auto expected_lower3 = lower3;
            auto test_lower3 = make_test_range(workspace.lower[3]);
            CHECK_EQUAL_RANGES(test_lower3, expected_lower3);
        }}
    };

    test(2,2,2,2);
    test(2,1,2,2);
    test(1,2,2,2);
    test(1,1,2,2);
}

BOOST_AUTO_TEST_CASE(test_peel_tomother_reverse) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1, int ploidy2, int ploidy3) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        const int sz2 = (ploidy2 == 1) ? 4 : 10;
        const int sz3 = (ploidy3 == 1) ? 4 : 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", ploidy2=" << ploidy2
                        << ", ploidy3=" << ploidy3
                        << ", test_num=" << test_num) {
            auto da = Mk::matrix(4, 1e-6, log(4));
            auto ma = Mk::matrix(4, 0.7e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), lower2(sz2), lower3(sz3), 
                                upper0(sz0), upper1(sz1),
                                expected_lower0(sz0), expected_lower1(sz1),
                                expected_super2(sz0*sz1), expected_upper2(sz2),
                                expected_super3(sz0*sz1), expected_upper3(sz3);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(lower2, [&](){ return xrand.get_double52(); });
            generate(lower3, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            generate(upper1, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(4);
            mats[2] = meiosis_matrix(ploidy0,da,ploidy1,ma);
            mats[3] = meiosis_matrix(ploidy0,ma,ploidy1,da);

            std::vector<double> p(sz0*sz1);

            // Calculate expected value
            for(int i=0;i<sz0;++i) {
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    p[ij] = 1.0;
                    double d = 0.0;
                    for(int k=0;k<sz2;++k) {           
                        d += mats[2](ij,k) * lower2[k];
                    }
                    p[ij] *= d;
                    d = 0.0;
                    for(int k=0;k<sz3;++k) {           
                        d += mats[3](ij,k) * lower3[k];
                    }
                    p[ij] *= d;
                }
            }

            for(int i=0;i<sz0;++i) {
                expected_lower0[i] = 0.0;
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    expected_lower0[i] += p[ij]*(upper1[j]*lower1[j]);
                }
                expected_lower0[i] *= lower0[i];
            }

            for(int i=0;i<sz0;++i) {
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    double d = 0.0;
                    for(int k=0;k<sz3;++k) {
                        d += mats[3](ij,k) * lower3[k];
                    }
                    expected_super2[ij] = d
                        *(upper0[i]*lower0[i])
                        *(lower1[j]*upper1[j]);
                }
            }

            for(int i=0;i<sz0;++i) {
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    double d = 0.0;
                    for(int k=0;k<sz2;++k) {           
                        d += mats[2](ij,k) * lower2[k];
                    }
                    expected_super3[ij] = d
                        *(upper0[i]*lower0[i])
                        *(lower1[j]*upper1[j]);
                }
            }

            for(int k=0;k<sz2;++k) {
                expected_upper2[k] = 0.0;
                for(int ij=0;ij<sz0*sz1;++ij) {
                    expected_upper2[k] += expected_super2[ij]*mats[2](ij,k);
                }
            }

            for(int k=0;k<sz3;++k) {
                expected_upper3[k] = 0.0;
                for(int ij=0;ij<sz0*sz1;++ij) {
                    expected_upper3[k] += expected_super3[ij]*mats[3](ij,k);
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(4);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.lower[2].resize(sz2);
            workspace.lower[3].resize(sz3);
            workspace.upper.resize(4);
            workspace.upper[0].resize(sz0);
            workspace.upper[1].resize(sz1);
            workspace.super.resize(4);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(lower2, workspace.lower[2].data());
            copy(lower3, workspace.lower[3].data());
            copy(upper0, workspace.upper[0].data());
            copy(upper1, workspace.upper[1].data());

            // do the peeling
            to_mother(workspace, {0,1,2,3}, mats);
            copy_n(workspace.lower[1].data(), sz1, expected_lower1.begin());
            to_mother_reverse(workspace, {0,1,2,3}, mats);

            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_CLOSE_RANGES(test_lower0, expected_lower0, prec);

            auto test_super2 = make_test_range(workspace.super[2]);
            CHECK_CLOSE_RANGES(test_super2, expected_super2, prec);

            auto test_super3 = make_test_range(workspace.super[3]);
            CHECK_CLOSE_RANGES(test_super3, expected_super3, prec);

            auto test_upper2 = make_test_range(workspace.upper[2]);
            CHECK_CLOSE_RANGES(test_upper2, expected_upper2, prec);

            auto test_upper3 = make_test_range(workspace.upper[3]);
            CHECK_CLOSE_RANGES(test_upper3, expected_upper3, prec);

            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_EQUAL_RANGES(test_lower1, expected_lower1);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto expected_upper1 = upper1;
            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_EQUAL_RANGES(test_upper1, expected_upper1);

            auto expected_lower2 = lower2;
            auto test_lower2 = make_test_range(workspace.lower[2]);
            CHECK_EQUAL_RANGES(test_lower2, expected_lower2);

            auto expected_lower3 = lower3;
            auto test_lower3 = make_test_range(workspace.lower[3]);
            CHECK_EQUAL_RANGES(test_lower3, expected_lower3);
        }}
    };

    test(2,2,2,2);
    test(2,1,2,2);
    test(1,2,2,2);
    test(1,1,2,2);
}

BOOST_AUTO_TEST_CASE(test_peel_tochild_reverse) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    auto test = [&xrand,prec](int ploidy0, int ploidy1, int ploidy2, int ploidy3) {
        const int sz0 = (ploidy0 == 1) ? 4 : 10;
        const int sz1 = (ploidy1 == 1) ? 4 : 10;
        const int sz2 = (ploidy2 == 1) ? 4 : 10;
        const int sz3 = (ploidy3 == 1) ? 4 : 10;
        for(int test_num=0; test_num < NUM_TEST; ++test_num) {
        BOOST_TEST_CONTEXT("ploidy0=" << ploidy0
                        << ", ploidy1=" << ploidy1
                        << ", ploidy2=" << ploidy2
                        << ", ploidy3=" << ploidy3
                        << ", test_num=" << test_num) {
            auto da = Mk::matrix(4, 1e-6, log(4));
            auto ma = Mk::matrix(4, 0.7e-6, log(4));

            // setup data
            std::vector<double> lower0(sz0), lower1(sz1), lower2(sz2), lower3(sz3), 
                                upper0(sz0), upper1(sz1),
                                expected_lower0(sz0), expected_lower1(sz1),
                                expected_super2(sz0*sz1), expected_upper2(sz2),
                                expected_super3(sz0*sz1), expected_upper3(sz3);
            generate(lower0, [&](){ return xrand.get_double52(); });
            generate(lower1, [&](){ return xrand.get_double52(); });
            generate(lower2, [&](){ return xrand.get_double52(); });
            generate(lower3, [&](){ return xrand.get_double52(); });
            generate(upper0, [&](){ return xrand.get_double52(); });
            generate(upper1, [&](){ return xrand.get_double52(); });
            TransitionMatrixVector mats(4);
            mats[2] = meiosis_matrix(ploidy0,da,ploidy1,ma);
            mats[3] = meiosis_matrix(ploidy0,ma,ploidy1,da);

            std::vector<double> p(sz0*sz1);

            // Calculate expected value
            for(int i=0;i<sz0;++i) {
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    p[ij] = 1.0;
                    double d = 0.0;
                    for(int k=0;k<sz2;++k) {           
                        d += mats[2](ij,k) * lower2[k];
                    }
                    p[ij] *= d;
                    d = 0.0;
                    for(int k=0;k<sz3;++k) {           
                        d += mats[3](ij,k) * lower3[k];
                    }
                    p[ij] *= d;
                }
            }

            for(int i=0;i<sz0;++i) {
                expected_lower0[i] = 0.0;
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    expected_lower0[i] += p[ij]*(upper1[j]*lower1[j]);
                }
                expected_lower0[i] *= lower0[i];
            }

            for(int j=0;j<sz1;++j) {
                expected_lower1[j] = 0.0;
                for(int i=0;i<sz0;++i) {
                    int ij = i*sz1+j;
                    expected_lower1[j] += p[ij]*(upper0[i]*lower0[i]);
                }
                expected_lower1[j] *= lower1[j];
            }

            for(int i=0;i<sz0;++i) {
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    double d = 0.0;
                    for(int k=0;k<sz3;++k) {
                        d += mats[3](ij,k) * lower3[k];
                    }
                    expected_super2[ij] = d
                        *(upper0[i]*lower0[i])
                        *(lower1[j]*upper1[j]);
                }
            }

            for(int i=0;i<sz0;++i) {
                for(int j=0;j<sz1;++j) {
                    int ij = i*sz1+j;
                    double d = 0.0;
                    for(int k=0;k<sz2;++k) {           
                        d += mats[2](ij,k) * lower2[k];
                    }
                    expected_super3[ij] = d
                        *(upper0[i]*lower0[i])
                        *(lower1[j]*upper1[j]);
                }
            }

            for(int k=0;k<sz2;++k) {
                expected_upper2[k] = 0.0;
                for(int ij=0;ij<sz0*sz1;++ij) {
                    expected_upper2[k] += expected_super2[ij]*mats[2](ij,k);
                }
            }

            for(int k=0;k<sz3;++k) {
                expected_upper3[k] = 0.0;
                for(int ij=0;ij<sz0*sz1;++ij) {
                    expected_upper3[k] += expected_super3[ij]*mats[3](ij,k);
                }
            }

            // setup the workspace and peeling operations
            workspace_t workspace;
            workspace.lower.resize(4);
            workspace.lower[0].resize(sz0);
            workspace.lower[1].resize(sz1);
            workspace.lower[2].resize(sz2);
            workspace.lower[3].resize(sz3);
            workspace.upper.resize(4);
            workspace.upper[0].resize(sz0);
            workspace.upper[1].resize(sz1);
            workspace.super.resize(4);
            copy(lower0, workspace.lower[0].data());
            copy(lower1, workspace.lower[1].data());
            copy(lower2, workspace.lower[2].data());
            copy(lower3, workspace.lower[3].data());
            copy(upper0, workspace.upper[0].data());
            copy(upper1, workspace.upper[1].data());

            // do the peeling
            to_child(workspace, {0,1,2,3}, mats);
            to_child_reverse(workspace, {0,1,2,3}, mats);

            auto test_lower0 = make_test_range(workspace.lower[0]);
            CHECK_CLOSE_RANGES(test_lower0, expected_lower0, prec);

            auto test_lower1 = make_test_range(workspace.lower[1]);
            CHECK_CLOSE_RANGES(test_lower1, expected_lower1, prec);

            auto test_super2 = make_test_range(workspace.super[2]);
            CHECK_CLOSE_RANGES(test_super2, expected_super2, prec);

            auto test_super3 = make_test_range(workspace.super[3]);
            CHECK_CLOSE_RANGES(test_super3, expected_super3, prec);

            auto test_upper2 = make_test_range(workspace.upper[2]);
            CHECK_CLOSE_RANGES(test_upper2, expected_upper2, prec);

            auto test_upper3 = make_test_range(workspace.upper[3]);
            CHECK_CLOSE_RANGES(test_upper3, expected_upper3, prec);

            auto expected_upper0 = upper0;
            auto test_upper0 = make_test_range(workspace.upper[0]);
            CHECK_EQUAL_RANGES(test_upper0, expected_upper0);

            auto expected_upper1 = upper1;
            auto test_upper1 = make_test_range(workspace.upper[1]);
            CHECK_EQUAL_RANGES(test_upper1, expected_upper1);

            auto expected_lower2 = lower2;
            auto test_lower2 = make_test_range(workspace.lower[2]);
            CHECK_EQUAL_RANGES(test_lower2, expected_lower2);

            auto expected_lower3 = lower3;
            auto test_lower3 = make_test_range(workspace.lower[3]);
            CHECK_EQUAL_RANGES(test_lower3, expected_lower3);
        }}
    };

    test(2,2,2,2);
    test(2,1,2,2);
    test(1,2,2,2);
    test(1,1,2,2);
}
