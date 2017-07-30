/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2016 Reed A. Cartwright
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

#include "../boost_test_helper.h"
#include "fixture_random_family.h"

using namespace dng;
using namespace dng::peel;
using dng::detail::make_test_range;

const int NUM_TEST = 100;
int g_seed_counter = 0;

BOOST_AUTO_TEST_CASE(test_peel_up_fast) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto m = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});

        // setup data
        std::vector<double> lower0(10), lower1(10), expected_lower0(10);
        generate(lower0, [&](){ return xrand.get_double52(); });
        generate(lower1, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(2);
        mats[1] = mitosis_diploid_matrix(m);

        // Calculated expected value
        for(int i=0;i<10;++i) {
            expected_lower0[i] = 0.0;
            for(int j=0;j<10;++j) {
                expected_lower0[i] += mats[1](i,j) * lower1[j];
            }
        }

        // setup the workspace and peeling operations
        workspace_t workspace;
        workspace.lower.resize(2);
        workspace.lower[0].resize(10);
        workspace.lower[1].resize(10);
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
}

BOOST_AUTO_TEST_CASE(test_peel_up) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto m = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});

        // setup data
        std::vector<double> lower0(10), lower1(10), expected_lower0(10);
        generate(lower0, [&](){ return xrand.get_double52(); });
        generate(lower1, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(2);
        mats[1] = mitosis_diploid_matrix(m);

        // Calculated expected value
        for(int i=0;i<10;++i) {
            expected_lower0[i] = 0.0;
            for(int j=0;j<10;++j) {
                expected_lower0[i] += mats[1](i,j) * lower1[j];
            }
            expected_lower0[i] *= lower0[i];
        }

        // setup the workspace and peeling operations
        workspace_t workspace;
        workspace.lower.resize(2);
        workspace.lower[0].resize(10);
        workspace.lower[1].resize(10);
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
}

BOOST_AUTO_TEST_CASE(test_peel_down_fast) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto m = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});

        // setup data
        std::vector<double> upper0(10), upper1(10), expected_upper1(10);
        generate(upper0, [&](){ return xrand.get_double52(); });
        generate(upper1, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(2);
        mats[1] = mitosis_diploid_matrix(m);

        // Calculated expected value
        for(int i=0;i<10;++i) {
            expected_upper1[i] = 0.0;
            for(int j=0;j<10;++j) {
                expected_upper1[i] += mats[1](j,i) * upper0[j];
            }
        }

        // setup the workspace and peeling operations
        workspace_t workspace;
        workspace.upper.resize(2);
        workspace.upper[0].resize(10);
        workspace.upper[1].resize(10);
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
}

BOOST_AUTO_TEST_CASE(test_peel_down) {
    const double prec = 4.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto m = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});

        // setup data
        std::vector<double> lower0(10), upper0(10), upper1(10), expected_upper1(10);
        generate(lower0, [&](){ return xrand.get_double52(); });
        generate(upper0, [&](){ return xrand.get_double52(); });
        generate(upper1, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(2);
        mats[1] = mitosis_diploid_matrix(m);

        // Calculated expected value
        for(int i=0;i<10;++i) {
            expected_upper1[i] = 0.0;
            for(int j=0;j<10;++j) {
                expected_upper1[i] += mats[1](j,i) * upper0[j]*lower0[j];
            }
        }

        // setup the workspace and peeling operations
        workspace_t workspace;
        workspace.lower.resize(1);
        workspace.lower[0].resize(10);
        copy(lower0, workspace.lower[0].data());
        workspace.upper.resize(2);
        workspace.upper[0].resize(10);
        workspace.upper[1].resize(10);
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
}

BOOST_AUTO_TEST_CASE(test_peel_tofather_fast) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto da = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});
        auto ma = f81::matrix(0.7e-6, {0.3, 0.1, 0.2, 0.4});

        // setup data
        std::vector<double> lower0(10), lower1(10), lower2(10), 
                upper1(10), expected_lower0(10);
        generate(lower0, [&](){ return xrand.get_double52(); });
        generate(lower1, [&](){ return xrand.get_double52(); });
        generate(lower2, [&](){ return xrand.get_double52(); });
        generate(upper1, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(3);
        mats[2] = meiosis_matrix(2,da,2,ma);

        // Calculated expected value
        for(int i=0;i<10;++i) {
            expected_lower0[i] = 0.0;
            for(int j=0;j<10;++j) {
                int ij = i*10+j;
                for(int k=0;k<10;++k) {
                    expected_lower0[i] += mats[2](ij,k) * lower2[k] * (lower1[j]*upper1[j]);
                }
            }
        }

        // setup the workspace and peeling operations
        workspace_t workspace;
        workspace.lower.resize(3);
        workspace.lower[0].resize(10);
        workspace.lower[1].resize(10);
        workspace.lower[2].resize(10);
        workspace.upper.resize(3);
        workspace.upper[1].resize(10);
        copy(lower0, workspace.lower[0].data());
        copy(lower1, workspace.lower[1].data());
        copy(lower2, workspace.lower[2].data());
        copy(upper1, workspace.upper[1].data());

        // do the peeling
        to_father_fast(workspace, {0,1,2}, mats);

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
    }}
}

BOOST_AUTO_TEST_CASE(test_peel_tofather) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto da = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});
        auto ma = f81::matrix(0.7e-6, {0.3, 0.1, 0.2, 0.4});

        // setup data
        std::vector<double> lower0(10), lower1(10), lower2(10), 
                upper1(10), expected_lower0(10);
        generate(lower0, [&](){ return xrand.get_double52(); });
        generate(lower1, [&](){ return xrand.get_double52(); });
        generate(lower2, [&](){ return xrand.get_double52(); });
        generate(upper1, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(3);
        mats[2] = meiosis_matrix(2,da,2,ma);

        // Calculated expected value
        for(int i=0;i<10;++i) {
            expected_lower0[i] = 0.0;
            for(int j=0;j<10;++j) {
                int ij = i*10+j;
                for(int k=0;k<10;++k) {
                    expected_lower0[i] += mats[2](ij,k) * lower2[k] * (lower1[j]*upper1[j]);
                }
            }
            expected_lower0[i] *= lower0[i];
        }

        // setup the workspace and peeling operations
        workspace_t workspace;
        workspace.lower.resize(3);
        workspace.lower[0].resize(10);
        workspace.lower[1].resize(10);
        workspace.lower[2].resize(10);
        workspace.upper.resize(3);
        workspace.upper[1].resize(10);
        copy(lower0, workspace.lower[0].data());
        copy(lower1, workspace.lower[1].data());
        copy(lower2, workspace.lower[2].data());
        copy(upper1, workspace.upper[1].data());

        // do the peeling
        to_father(workspace, {0,1,2}, mats);

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
    }}
}

BOOST_AUTO_TEST_CASE(test_peel_tomother_fast) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto da = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});
        auto ma = f81::matrix(0.7e-6, {0.3, 0.1, 0.2, 0.4});

        // setup data
        std::vector<double> lower0(10), lower1(10), lower2(10), 
                upper0(10), expected_lower1(10);
        generate(lower0, [&](){ return xrand.get_double52(); });
        generate(lower1, [&](){ return xrand.get_double52(); });
        generate(lower2, [&](){ return xrand.get_double52(); });
        generate(upper0, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(3);
        mats[2] = meiosis_matrix(2,da,2,ma);

        // Calculated expected value
        for(int i=0;i<10;++i) {
            expected_lower1[i] = 0.0;
            for(int j=0;j<10;++j) {
                int ij = j*10+i;
                for(int k=0;k<10;++k) {
                    expected_lower1[i] += mats[2](ij,k) * lower2[k] * (lower0[j]*upper0[j]);
                }
            }
        }

        // setup the workspace and peeling operations
        workspace_t workspace;
        workspace.lower.resize(3);
        workspace.lower[0].resize(10);
        workspace.lower[1].resize(10);
        workspace.lower[2].resize(10);
        workspace.upper.resize(3);
        workspace.upper[0].resize(10);
        copy(lower0, workspace.lower[0].data());
        copy(lower1, workspace.lower[1].data());
        copy(lower2, workspace.lower[2].data());
        copy(upper0, workspace.upper[0].data());

        // do the peeling
        to_mother_fast(workspace, {0,1,2}, mats);

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
    }}
}

BOOST_AUTO_TEST_CASE(test_peel_tomother) {
    const double prec = 8.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto da = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});
        auto ma = f81::matrix(0.7e-6, {0.3, 0.1, 0.2, 0.4});

        // setup data
        std::vector<double> lower0(10), lower1(10), lower2(10), 
                upper0(10), expected_lower1(10);
        generate(lower0, [&](){ return xrand.get_double52(); });
        generate(lower1, [&](){ return xrand.get_double52(); });
        generate(lower2, [&](){ return xrand.get_double52(); });
        generate(upper0, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(3);
        mats[2] = meiosis_matrix(2,da,2,ma);

        // Calculated expected value
        for(int i=0;i<10;++i) {
            expected_lower1[i] = 0.0;
            for(int j=0;j<10;++j) {
                int ij = j*10+i;
                for(int k=0;k<10;++k) {
                    expected_lower1[i] += mats[2](ij,k) * lower2[k] * (lower0[j]*upper0[j]);
                }
            }
            expected_lower1[i] *= lower1[i];
        }

        // setup the workspace and peeling operations
        workspace_t workspace;
        workspace.lower.resize(3);
        workspace.lower[0].resize(10);
        workspace.lower[1].resize(10);
        workspace.lower[2].resize(10);
        workspace.upper.resize(3);
        workspace.upper[0].resize(10);
        copy(lower0, workspace.lower[0].data());
        copy(lower1, workspace.lower[1].data());
        copy(lower2, workspace.lower[2].data());
        copy(upper0, workspace.upper[0].data());

        // do the peeling
        to_mother(workspace, {0,1,2}, mats);

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
    }}
}

BOOST_AUTO_TEST_CASE(test_peel_tochild_fast) {
    const double prec = 16.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto da = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});
        auto ma = f81::matrix(0.7e-6, {0.3, 0.1, 0.2, 0.4});

        // setup data
        std::vector<double> lower0(10), lower1(10), lower2(10), 
                            upper0(10), upper1(10), upper2(10),
                            expected_upper2(10);
        generate(lower0, [&](){ return xrand.get_double52(); });
        generate(lower1, [&](){ return xrand.get_double52(); });
        generate(lower2, [&](){ return xrand.get_double52(); });
        generate(upper0, [&](){ return xrand.get_double52(); });
        generate(upper1, [&](){ return xrand.get_double52(); });
        generate(upper2, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(3);
        mats[2] = meiosis_matrix(2,da,2,ma);

        // Calculated expected value
        for(int k=0;k<10;++k) {
            expected_upper2[k] = 0.0;
            for(int i=0,ij=0;i<10;++i) {
                for(int j=0;j<10;++j,++ij) {
                    expected_upper2[k] += mats[2](ij,k) * (lower0[i]*upper0[i]) * (lower1[j]*upper1[j]);
                }
            }
        }

        // setup the workspace and peeling operations
        workspace_t workspace;
        workspace.lower.resize(3);
        workspace.lower[0].resize(10);
        workspace.lower[1].resize(10);
        workspace.lower[2].resize(10);
        workspace.upper.resize(3);
        workspace.upper[0].resize(10);
        workspace.upper[1].resize(10);
        workspace.upper[2].resize(10);
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
}

BOOST_AUTO_TEST_CASE(test_peel_tochild) {
    const double prec = 32.0*DBL_EPSILON;
    using boost::copy;
    using boost::generate;

    xorshift64 xrand(++g_seed_counter);

    for(int test_num=0; test_num < NUM_TEST; ++test_num) {
    BOOST_TEST_CONTEXT("test_num=" << test_num) {
        auto da = f81::matrix(1e-6, {0.4, 0.2, 0.1, 0.3});
        auto ma = f81::matrix(0.7e-6, {0.3, 0.1, 0.2, 0.4});

        // setup data
        std::vector<double> lower0(10), lower1(10), lower2(10), lower3(10), 
                            upper0(10), upper1(10), upper2(10), upper3(10),
                            expected_upper2(10);
        generate(lower0, [&](){ return xrand.get_double52(); });
        generate(lower1, [&](){ return xrand.get_double52(); });
        generate(lower2, [&](){ return xrand.get_double52(); });
        generate(lower3, [&](){ return xrand.get_double52(); });
        generate(upper0, [&](){ return xrand.get_double52(); });
        generate(upper1, [&](){ return xrand.get_double52(); });
        generate(upper2, [&](){ return xrand.get_double52(); });
        generate(upper3, [&](){ return xrand.get_double52(); });
        TransitionMatrixVector mats(4);
        mats[2] = meiosis_matrix(2,da,2,ma);
        mats[3] = meiosis_matrix(2,da,2,ma);

        // Calculated expected value
        for(int k=0;k<10;++k) {
            expected_upper2[k] = 0.0;
            for(int i=0,ij=0;i<10;++i) {
                for(int j=0;j<10;++j,++ij) {
                    for(int h=0;h<10;++h) {
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
        workspace.lower[0].resize(10);
        workspace.lower[1].resize(10);
        workspace.lower[2].resize(10);
        workspace.lower[3].resize(10);
        workspace.upper.resize(4);
        workspace.upper[0].resize(10);
        workspace.upper[1].resize(10);
        workspace.upper[2].resize(10);
        workspace.upper[3].resize(10);
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
}
