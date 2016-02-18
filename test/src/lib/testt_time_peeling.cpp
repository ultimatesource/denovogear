//
// Created by steven on 1/15/16.
//

#define BOOST_TEST_MODULE dng::lib::peeling

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/timer/timer.hpp>

#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <dng/peeling.h>
#include <chrono>

#include "boost_test_helper.h"
#include "testh_helper_peeling.h"

//Time trail/perforamnce test. Will take longer to run than normal test.
using namespace dng;
namespace utf = boost::unit_test;

void to_father_fast_original(dng::peel::workspace_t &work, const dng::peel::family_members_t &family,
                                        const TransitionVector &mat) {

    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Include Mom
    work.paired_buffer.resize(10, 10);
    work.lower[dad] = (work.paired_buffer.matrix() * (work.upper[mom] *
                                                      work.lower[mom]).matrix()).array();
    work.paired_buffer.resize(100, 1);
}


void to_father_original(dng::peel::workspace_t &work, const dng::peel::family_members_t &family,
                                   const TransitionVector &mat) {

    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Include Mom
    work.paired_buffer.resize(10, 10);
    work.lower[dad] *= (work.paired_buffer.matrix() * (work.upper[mom] *
                                                       work.lower[mom]).matrix()).array();
    work.paired_buffer.resize(100, 1);
}


BOOST_FIXTURE_TEST_SUITE(test_peeling_speed, Fx)

BOOST_AUTO_TEST_CASE(test_speed_to_father) {

    const int num_time_trial = 1e3;
    const int num_repeat = 1e4;

    std::chrono::duration<double> elapsed_seconds, elapsed_update, elapsed_original;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::time_t start_time, end_time;

    elapsed_update = {};
    elapsed_original = {};
    int verbose_gap = num_time_trial/10;
    for (int r = 0; r < num_time_trial; ++r) {
        if( (r % verbose_gap) == 0){
            std::cout << "At "<< (r/verbose_gap) << "0 %"<< std::endl;
        }


        init_family();

        dng::peel::workspace_t workspace;
        dng::TransitionVector full_matrix;
        copy_family_to_workspace(workspace, full_matrix, total_family_size,
                                     lower_array, upper_array, trans_matrix);

        GenotypeArray temp_result;

        {
            //Updated version
//            boost::timer::auto_cpu_timer measure_speed(std::cerr);
            temp_result = GenotypeArray::Zero();
            start = std::chrono::system_clock::now();
            for (int t = 0; t < num_repeat; ++t) {
                workspace.lower[0] = lower_array[0];
                //            dng::peel::to_father_fast(workspace, family, full_matrix);
                dng::peel::to_father(workspace, family, full_matrix);
                GenotypeArray result = workspace.lower[0];
                temp_result += result;
            }
            end = std::chrono::system_clock::now();
            elapsed_seconds = end - start;
            elapsed_update += elapsed_seconds;
            //        std::cout << "Updated time: " << elapsed_seconds.count() << "s\t" << temp_result.sum() << std::endl;
        }

        {
            //Original version
//            boost::timer::auto_cpu_timer measure_speed(std::cerr);
            temp_result = GenotypeArray::Zero();
            start = std::chrono::system_clock::now();
            for (int t = 0; t < num_repeat; ++t) {
                workspace.lower[0] = lower_array[0];
//            dng::peel::to_father_fast_original(workspace, family, full_matrix);
                to_father_original(workspace, family, full_matrix);
                GenotypeArray result_original = workspace.lower[0];
                temp_result += result_original;
            }
            end = std::chrono::system_clock::now();
            elapsed_seconds = end - start;
            elapsed_original += elapsed_seconds;
//        std::cout << "Original time: " << elapsed_seconds.count() << "s\t" << temp_result.sum() << std::endl;
        }
    }

    double total = num_time_trial*num_repeat;
    std::cout << "Summary: "<< total << " calls." << std::endl;
    std::cout << "Original total: " << elapsed_original.count() <<
            " each: " << elapsed_original.count()/total << std::endl;
    std::cout << "Update   total: " << elapsed_update.count() <<
            " each: " << elapsed_update.count()/total << std::endl;
    std::cout << "delta/origin: "<< (elapsed_update-elapsed_original)/elapsed_original << std::endl;
}


BOOST_AUTO_TEST_SUITE_END()



/*
Summary: 1e+07 calls.
*** No errors detected
Original total: 34.8541 each: 3.48541e-06
Update   total: 33.3654 each: 3.33654e-06
delta/origin: -0.0427128

*/
