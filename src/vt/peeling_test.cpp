//
// Created by steven on 1/15/16.
//

#define BOOST_TEST_MODULE dng::lib::peeling

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <string>
#include <cstdlib>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <ctime>

#include <dng/peeling.h>
#include "assert_helper.h"

//#include <boost/test/data/test_case.hpp>
////#include <boost/test/data/monomorphic.hpp>
//namespace data = boost::unit_test::data;

using namespace dng;
namespace utf = boost::unit_test;


int num_test = 10;

std::random_device rd;
std::mt19937 gen(rd());

struct Fx {
    std::string s;

    const int child_offset = 2;


    dng::peel::family_members_t family;

    std::vector<TransitionMatrix> m;
    std::vector<GenotypeArray> u;
    std::vector<GenotypeArray> g;

    int num_child;
    int total_family_size;
    std::uniform_int_distribution<> rand_unif;
//    std::uniform_int_distribution<> rand_unif(1, 10);
//    std::uniform_int_distribution<> rand_unif(1, 10);

    Fx(std::string s = "") : s(s) {
//        std::srand(std::time(0));
//        std::mt19937 gen(rd2());

        BOOST_TEST_MESSAGE("set up " << s);

        rand_unif = std::uniform_int_distribution<>(1,10);
        init_family();
    }

    void init_family(){



        int num_child = rand_unif(gen);

        total_family_size = num_child + child_offset;
        std::cout << num_child << std::endl;

        family.clear();
        m.resize(total_family_size);
        u.resize(total_family_size);
        g.resize(total_family_size);

        for (int k = 0; k < total_family_size; ++k) {
            family.push_back(k);
            g[k] = GenotypeArray::Random();
            if (k < child_offset) {
                m[k] = TransitionMatrix::Random(10, 10);
                u[k] = GenotypeArray::Random();
            }
            else {
                m[k] = TransitionMatrix::Random(100, 10);
            }
        }
    }

    ~Fx() {
        BOOST_TEST_MESSAGE("tear down " << s);
    }
};

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }


// TODO: Example of BOOST_DATA_TEST_CASE and BOOST_PARAM_TEST_CASE.
// TODO: Should be able to replace the for loop with these.
// TODO: Might not be able to use fixture.
//
//
//std::vector<int> test_types;//
//
//BOOST_AUTO_TEST_SUITE(suite1,
//  * utf::fixture<Fx>(std::string("FX"))
//  * utf::fixture<Fx>(std::string("FX2")))
//
//  BOOST_AUTO_TEST_CASE(test1, * utf::fixture(&setup, &teardown))
//  {
//    BOOST_TEST_MESSAGE("running test1");
//    BOOST_TEST(true);
//  }
//
//  BOOST_AUTO_TEST_CASE(test2)
//  {
//    BOOST_TEST_MESSAGE("running test2");
//    BOOST_TEST(true);
//  }
//
//BOOST_DATA_TEST_CASE( test_case_arity1, data::xrange(5), my_var )
//{
//    BOOST_TEST_MESSAGE("running data: ");
//    BOOST_TEST((my_var <= 4 && my_var >= 0));
//}
//BOOST_PARAM_TEST_CASE(test_function, params_begin, params_end);
//BOOST_AUTO_TEST_SUITE_END()
//


BOOST_FIXTURE_TEST_SUITE(test_peeling_suite, Fx)
//BOOST_AUTO_TEST_SUITE(test_peeling_suite,  * utf::fixture<Fx>(std::string("FX")) )

    BOOST_AUTO_TEST_CASE(test_sum_over_child, * utf::fixture(&setup, &teardown)) {

        for (int t = 0; t < num_test; ++t) {

//            peel::family_members_t family;
            init_family();
//            std::uniform_int_distribution<> rand_unif(1, 10);
//            int num_child = rand_unif(gen);
//            int total_family_size = num_child + child_offset;
//
//            std::vector<TransitionMatrix> m(total_family_size);
////        std::vector<GenotypeArray> u(total_family_size);
//            std::vector<GenotypeArray> g(total_family_size);
//
//            for (int k = 0; k < total_family_size; ++k) {
//                family.push_back(k);
//                g[k] = GenotypeArray::Random();
//                if (k < child_offset) {
//                    m[k] = TransitionMatrix::Random(10, 10);
////                u[k] = GenotypeArray::Random();
//                }
//                else {
//                    m[k] = TransitionMatrix::Random(100, 10);
//                }
//            }

            PairedGenotypeArray expected = PairedGenotypeArray::Ones(100, 1);
            for (int k = 2; k < total_family_size; ++k) {
                PairedGenotypeArray temp_array = PairedGenotypeArray::Zero(100, 1);

                for (int i = 0; i < m[k].rows(); ++i) {
                    double x = 0;
                    for (int j = 0; j < m[k].cols(); ++j) {
                        x += m[k](i, j) * g[k][j];
                    }
                    temp_array(i, 0) = x;
                }
                expected *= temp_array;
            }

            dng::peel::workspace_t workspace;
            workspace.Resize(total_family_size);
            dng::TransitionVector full_matrix(total_family_size);
            for (int l = 0; l < total_family_size; ++l) {
                workspace.lower[l] = g[l];
//            workspace.upper[l] = u[l];
                full_matrix[l] = m[l];
            }


            PairedGenotypeArray result = dng::peel::sum_over_child(workspace, family, full_matrix);
            for (int i = 0; i < 100; ++i) {
                BOOST_CHECK_CLOSE(expected(i, 0), result(i, 0), BOOST_CLOSE_THRESHOLD);
            }
        }


    }


    BOOST_AUTO_TEST_CASE(test_up_core) {

        peel::family_members_t family1{0, 1}; //0 parent, 1 child
        dng::peel::workspace_t workspace;

        for (int t = 0; t < num_test; ++t) {
            const TransitionMatrix m = TransitionMatrix::Random(10, 10);
            const GenotypeArray g = GenotypeArray::Random();
            GenotypeArray expected = GenotypeArray::Zero();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    expected[i] += m(i, j) * g[j];
                }
            }

            workspace.Resize(2);
            dng::TransitionVector full_matrix{2};
            full_matrix[0] = {};
            full_matrix[1] = m;
            workspace.lower[1] = g;

            GenotypeArray result = dng::peel::up_core(workspace, family1, full_matrix);
            for (int i = 0; i < 10; ++i) {
                BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
            }
        }

    }

    BOOST_AUTO_TEST_CASE(test_up) {

        peel::family_members_t family1{0, 1}; //0 parent, 1 child
        dng::peel::workspace_t workspace;

        for (int t = 0; t < num_test; ++t) {
            const TransitionMatrix m = TransitionMatrix::Random(10, 10);
            const GenotypeArray g0 = GenotypeArray::Random();
            const GenotypeArray g1 = GenotypeArray::Random();
            GenotypeArray expected = GenotypeArray::Zero();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    expected[i] += m(i, j) * g1[j];
                }
                expected[i] *= g0[i];
            }

            workspace.Resize(2);
            dng::TransitionVector full_matrix{2};
            full_matrix[0] = {};
            full_matrix[1] = m;
            workspace.lower[0] = g0;
            workspace.lower[1] = g1;

            dng::peel::up(workspace, family1, full_matrix);
            GenotypeArray result = workspace.lower[0];
            for (int i = 0; i < 10; ++i) {
                BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
            }
        }

    }


    BOOST_AUTO_TEST_CASE(test_to_father) {

//        int child_offset = 2;
//        std::uniform_int_distribution<> rand_unif(1, 10);

        for (int t = 0; t < num_test; ++t) {

//            peel::family_members_t family;
//
//            int num_child = rand_unif(gen);
//
//
//            int total_family_size = num_child + child_offset;
//            std::vector<TransitionMatrix> m(total_family_size);
//            std::vector<GenotypeArray> u(total_family_size);
//            std::vector<GenotypeArray> g(total_family_size);
//
//            for (int k = 0; k < total_family_size; ++k) {
//                family.push_back(k);
//                g[k] = GenotypeArray::Random();
//                if (k < child_offset) {
//                    m[k] = TransitionMatrix::Random(10, 10);
//                    u[k] = GenotypeArray::Random();
//                }
//                else {
//                    m[k] = TransitionMatrix::Random(100, 10);
//                }
//            }

            init_family();

            PairedGenotypeArray all_child = PairedGenotypeArray::Ones(100, 1);
            for (int c = child_offset; c < total_family_size; ++c) {
                PairedGenotypeArray one_child = PairedGenotypeArray::Zero(100, 1);
                for (int i = 0; i < 100; ++i) {
                    for (int j = 0; j < 10; ++j) {
                        one_child(i, 0) += m[c](i, j) * g[c][j];
                    }
                }
                all_child *= one_child;
            }
            all_child.resize(10, 10);

            GenotypeArray ga_mum;
            for (int j = 0; j < 10; ++j) {
                ga_mum[j] = u[1][j] * g[1][j];
            }

            GenotypeArray expected = GenotypeArray::Zero();
            GenotypeArray expected_fast = GenotypeArray::Zero();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    expected[i] += all_child(i, j) * ga_mum[j];
                }
                expected_fast[i] = expected[i];
                expected[i] *= g[0][i];
            }
//        std::cout << expected.sum() << std::endl;


            dng::peel::workspace_t workspace;
            workspace.Resize(total_family_size);
            dng::TransitionVector full_matrix(total_family_size);
            for (int l = 0; l < total_family_size; ++l) {
                workspace.lower[l] = g[l];
                workspace.upper[l] = u[l];
                full_matrix[l] = m[l];
            }

            dng::peel::to_father(workspace, family, full_matrix);
            GenotypeArray result = workspace.lower[0];
            for (int i = 0; i < 10; ++i) {
                BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
            }

            dng::peel::to_father_fast(workspace, family, full_matrix);
            GenotypeArray result_fast = workspace.lower[0];
            for (int i = 0; i < 10; ++i) {
                BOOST_CHECK_CLOSE(expected_fast[i], result_fast[i], BOOST_CLOSE_THRESHOLD);
            }


        }

    }


    BOOST_AUTO_TEST_CASE(test_to_mother) {

//        int child_offset = 2;
//        std::uniform_int_distribution<> rand_unif(1, 10);


        for (int t = 0; t < num_test; ++t) {

//            int num_child = rand_unif(gen);
//            int total_family_size = num_child + child_offset;
//            std::vector<TransitionMatrix> m(total_family_size);
//            std::vector<GenotypeArray> u(total_family_size);
//            std::vector<GenotypeArray> g(total_family_size);
//
//            peel::family_members_t family;
//            for (int k = 0; k < total_family_size; ++k) {
//                family.push_back(k);
//                g[k] = GenotypeArray::Random();
//                if (k < child_offset) {
//                    m[k] = TransitionMatrix::Random(10, 10);
//                    u[k] = GenotypeArray::Random();
//                }
//                else {
//                    m[k] = TransitionMatrix::Random(100, 10);
//                }
//            }

            init_family();

            PairedGenotypeArray all_child = PairedGenotypeArray::Ones(100, 1);
            for (int c = child_offset; c < total_family_size; ++c) {
                PairedGenotypeArray one_child = PairedGenotypeArray::Zero(100, 1);
                for (int i = 0; i < 100; ++i) {
                    for (int j = 0; j < 10; ++j) {
                        one_child(i, 0) += m[c](i, j) * g[c][j];
                    }
                }
                all_child *= one_child;
            }
            all_child.resize(10, 10);

            GenotypeArray ga_dad;
            for (int j = 0; j < 10; ++j) {
                ga_dad[j] = u[0][j] * g[0][j];
            }

            GenotypeArray expected = GenotypeArray::Zero();
            GenotypeArray expected_fast = GenotypeArray::Zero();
            auto temp_child = all_child.transpose();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    expected[i] += temp_child(i, j) * ga_dad[j];
                }
                expected_fast[i] = expected[i];
                expected[i] *= g[1][i];
            }
//        std::cout << expected.sum() << std::endl;

            dng::peel::workspace_t workspace;
            workspace.Resize(total_family_size);
            dng::TransitionVector full_matrix(total_family_size);
            for (int l = 0; l < total_family_size; ++l) {
                workspace.lower[l] = g[l];
                workspace.upper[l] = u[l];
                full_matrix[l] = m[l];
            }


            dng::peel::to_mother(workspace, family, full_matrix);
            GenotypeArray result = workspace.lower[1];
            for (int i = 0; i < 10; ++i) {
                BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
            }

            dng::peel::to_mother_fast(workspace, family, full_matrix);
            GenotypeArray result_fast = workspace.lower[1];
            for (int i = 0; i < 10; ++i) {
                BOOST_CHECK_CLOSE(expected_fast[i], result_fast[i], BOOST_CLOSE_THRESHOLD);
            }


        }

    }


    BOOST_AUTO_TEST_CASE(test_to_child) {

//        int child_offset = 2;
//        std::uniform_int_distribution<> rand_unif(2, 10);

        rand_unif = std::uniform_int_distribution<>(2,10);

        for (int t = 0; t < num_test; ++t) {

//            peel::family_members_t family;
//
//            int num_child = rand_unif(gen);
//            int total_family_size = num_child + child_offset;
//            std::vector<TransitionMatrix> m(total_family_size);
//            std::vector<GenotypeArray> u(total_family_size);
//            std::vector<GenotypeArray> g(total_family_size);
//
//            for (int k = 0; k < total_family_size; ++k) {
//                family.push_back(k);
//                g[k] = GenotypeArray::Random();
//                if (k < child_offset) {
//                    m[k] = TransitionMatrix::Random(10, 10);
//                    u[k] = GenotypeArray::Random();
//                }
//                else {
//                    m[k] = TransitionMatrix::Random(100, 10);
//                }
//            }

            init_family();
            
            PairedGenotypeArray all_child = PairedGenotypeArray::Ones(100, 1);
            for (int c = child_offset + 1; c < total_family_size; ++c) {
                PairedGenotypeArray one_child = PairedGenotypeArray::Zero(100, 1);
                for (int i = 0; i < 100; ++i) {
                    for (int j = 0; j < 10; ++j) {
                        one_child(i, 0) += m[c](i, j) * g[c][j];
                    }
                }
                all_child *= one_child;
            }

            GenotypeArray dad_array;
            GenotypeArray mom_array;
            for (int j = 0; j < 10; ++j) {
                dad_array[j] = u[0][j] * g[0][j];
                mom_array[j] = u[1][j] * g[1][j];
            }
            PairedGenotypeArray parents_array = PairedGenotypeArray::Zero(100, 1);
            for (int d = 0; d < 10; ++d) {
                for (int m = 0; m < 10; ++m) {
                    parents_array(d * 10 + m, 0) = dad_array[d] * mom_array[m];
                }
            }
            PairedGenotypeArray parents_children = parents_array * all_child;

            GenotypeArray expected = GenotypeArray::Zero();
            GenotypeArray expected_fast = GenotypeArray::Zero();
            auto temp_m = m[2].transpose();
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 100; ++j) {
                    expected[i] += temp_m(i, j) * parents_children(j, 0);
                    expected_fast[i] += temp_m(i, j) * parents_array(j, 0);
                }
            }


            dng::peel::workspace_t workspace;
            workspace.Resize(total_family_size);
            dng::TransitionVector full_matrix(total_family_size);
            for (int l = 0; l < total_family_size; ++l) {
                workspace.lower[l] = g[l];
                workspace.upper[l] = u[l];
                full_matrix[l] = m[l];
            }


            dng::peel::to_child(workspace, family, full_matrix);
            GenotypeArray result = workspace.upper[2];
            for (int i = 0; i < expected.size(); ++i) {
                BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
            }

            peel::family_members_t family_fast{family[0], family[1], family[2]};
            dng::peel::to_child_fast(workspace, family_fast, full_matrix);
            GenotypeArray result_fast = workspace.upper[2];
            for (int i = 0; i < expected_fast.size(); ++i) {
                BOOST_CHECK_CLOSE(expected_fast[i], result_fast[i], BOOST_CLOSE_THRESHOLD);

            }


        }


    }

BOOST_AUTO_TEST_SUITE_END()
