#define BOOST_TEST_MODULE dng::genotype

#include <dng/genotyper.h>

#include "../testing.h"

using namespace dng;
using namespace dng::genotype;

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
