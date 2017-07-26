#define BOOST_TEST_MODULE dng::genotype

#include <dng/genotyper.h>

#include "../testing.h"

using namespace dng;
using namespace dng::genotype;
using dng::detail::make_test_range;
using d4 = std::array<double,4>;

BOOST_AUTO_TEST_CASE(test_diploid_nucleotides) {
    for(int a : {0,1,2,3}) {
        for(int b : {0,1,2,3}) {
            BOOST_TEST_CONTEXT("a=" << a << ", b=" << b) {
                int gt = folded_diploid_genotypes_matrix[a][b];
                
                BOOST_CHECK((folded_diploid_nucleotides[gt][0] == a && folded_diploid_nucleotides[gt][1] == b) ||
                    folded_diploid_nucleotides[gt][0] == b && folded_diploid_nucleotides[gt][1] == a);
                int unfolded = a*4+b;
                int gt2 = folded_diploid_genotypes[unfolded];
                BOOST_CHECK((folded_diploid_nucleotides[gt2][0] == a && folded_diploid_nucleotides[gt2][1] == b) ||
                    folded_diploid_nucleotides[gt2][0] == b && folded_diploid_nucleotides[gt2][1] == a);

                if(a <= b) {
                    BOOST_CHECK_EQUAL(unfolded_diploid_genotypes_upper[gt2], unfolded);
                }
                if(b <= a) {
                    BOOST_CHECK_EQUAL(unfolded_diploid_genotypes_lower[gt2], unfolded);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_log_pochhammer) {
    double prec = 32*DBL_EPSILON;

    auto test = [prec](double a) -> void {
        BOOST_TEST_CONTEXT("a=" << a)
    {
        using dng::genotype::detail::log_pochhammer;
        log_pochhammer pochy(a);
        std::vector<double> test_values, expected_values;
        double d = 0.0;
        for(int i=0; i < 1000; ++i) {
            d += log(i+a);
            expected_values.push_back(d);
            test_values.push_back(pochy(i+1));
        }
        CHECK_CLOSE_RANGES( test_values, expected_values, prec);
    }
    };

    test(1e0);
    test(1e2);
    test(1e4);
    test(1e6);
    test(1e8);

    test(1e-2);
    test(1e-4);
    test(1e-6);
    test(1e-8);

    double low_phi = DBL_EPSILON/2.0;
    test((1.0 - low_phi) / low_phi);
}