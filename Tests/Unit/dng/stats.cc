#define BOOST_TEST_MODULE dng::stats

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <array>
#include <numeric>
#include <limits>

#include <dng/stats.h>

using namespace std;

bool check_sum(double a, double b ) {
    if(a == b) {
        return true;
    }
    cerr << setprecision(std::numeric_limits<double>::max_digits10);
    cerr << "  ERROR: Values do not match. Expected '" << b << "'. Got '"
         << a << "'.\n";
    return false;
}

// use many tests borrowed from python
BOOST_AUTO_TEST_CASE(test_exact_sum) {
    using dng::stats::ExactSum;
    using dng::stats::exact_sum;
    BOOST_CHECK(ExactSum() == 0.0);
    {
        ExactSum sum(1.0);
        BOOST_CHECK(check_sum(sum,1.0));
        sum += 1.0;
        BOOST_CHECK(check_sum(sum,2.0));
        sum += -1.0;
        BOOST_CHECK(check_sum(sum,1.0));
        sum(1.0);
        BOOST_CHECK(check_sum(sum,2.0));
        sum(-1.0);
        BOOST_CHECK(check_sum(sum,1.0));
        (sum(1.0) += 1.0).add(-1.0) += -1.0;
        BOOST_CHECK(sum(-1.0).result() == 0.0);
        BOOST_CHECK(sum == false && sum.failed() == false);
    }
    BOOST_CHECK(check_sum(exact_sum({0.0}),0.0));
    BOOST_CHECK(check_sum(exact_sum(std::array<double,2>({1.0,-1.0})), 0.0));
    BOOST_CHECK(check_sum(exact_sum({1e100, 1.0, -1e100, 1e-100, 1e50, -1.0, -1e50}), 1e-100));
    BOOST_CHECK(check_sum(exact_sum({pow(2.0,53), -0.5, -pow(2.0,-54)}), pow(2.0,53)-1.0));
    BOOST_CHECK(check_sum(exact_sum({pow(2.0,53), 1.0, pow(2.0,-100)}), pow(2.0,53)+2.0));
    BOOST_CHECK(check_sum(exact_sum({pow(2.0,53)+10.0, 1.0, pow(2.0,-100)}), pow(2.0,53)+12.0));
    BOOST_CHECK(check_sum(exact_sum({pow(2.0,53)-4.0, 0.5, pow(2.0,-54)}), pow(2.0,53)-3.0));
    BOOST_CHECK(check_sum(exact_sum({1e16, 1.0, 1e-16}), 10000000000000002.0));
    BOOST_CHECK(check_sum(exact_sum({1e16-2.0, 1.0-pow(2.0,-53), -(1e16-2.0), -(1.0-pow(2.0,-53))}), 0.0));
    {
        ExactSum sum;
        for(int i=1;i<1001;++i) {
            sum(1.0/i);
        }
        BOOST_CHECK(check_sum(sum,atof("0x1.df11f45f4e61ap+2")));
    }
    {
        ExactSum sum;
        for(int i=1;i<1001;++i) {
            sum(pow(-1.0,i)/i);
        }
        BOOST_CHECK(check_sum(sum,atof("-0x1.62a2af1bd3624p-1")));
    }
    {
        ExactSum sum;
        for(int i=0;i<1000;++i) {
            sum(pow(1.7,i+1)-pow(1.7,i));
        }
        sum(-pow(1.7,1000));
        BOOST_CHECK(check_sum(sum,-1.0));
    }
    {
        ExactSum sum;
        for(int i=-1074; i < 972; i += 2) {
            sum(pow(2.0,i) - pow(2.0,i+50) + pow(2.0,i+52));
        }
        sum(-pow(2.0,1022));
        BOOST_CHECK(check_sum(sum,atof("0x1.5555555555555p+970")));
    }
    BOOST_CHECK(std::isnan(exact_sum({1.0,(double)NAN})));
    BOOST_CHECK(std::isnan(exact_sum({(double)-INFINITY,(double)INFINITY})));

    BOOST_CHECK(check_sum(exact_sum({pow(2,1023),pow(2,1023)}),HUGE_VAL));
    BOOST_CHECK(check_sum(exact_sum({1.0, (double)INFINITY}),INFINITY));
    BOOST_CHECK(check_sum(exact_sum({1.0, (double)-INFINITY}),-INFINITY));
}


// Really the exact test (no approx.), same value as returned by R
BOOST_AUTO_TEST_CASE(test_fisher_few_reads){
    BOOST_CHECK_CLOSE(dng::stats::fisher_exact_test(6,12,12,5),0.04371017, 0.00001);
}


// When ther are > 512 total reads a g-test is used to test approximate the
// Fisher Test. Test value calculated in R.
BOOST_AUTO_TEST_CASE(test_fisher_many_reads){
    BOOST_CHECK_CLOSE(dng::stats::fisher_exact_test(20, 200, 19, 301), 0.1679341, 0.0001);

}

//Test value for scipy k-sample test
BOOST_AUTO_TEST_CASE(test_AD){
    std::vector<int> a = {40, 31, 35, 40, 40, 32, 33};           
    std::vector<int> b = {21, 31, 33, 34, 34, 40, 42, 20} ;
    BOOST_CHECK_CLOSE(dng::stats::ad_two_sample_test(a,b), -0.52579911592960638, 0.00001);
}


