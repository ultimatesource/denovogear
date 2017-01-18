#define BOOST_TEST_MODULE dng::genotype

#include <dng/genotyper.h>

#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

namespace data = boost::unit_test::data;

using namespace dng;
using namespace dng::genotype;

BOOST_DATA_TEST_CASE(test_diploid_nucleotides, data::xrange(4) * data::xrange(4), a, b)
{
	int gt = folded_diploid_genotypes_matrix[a][b];
	BOOST_CHECK((folded_diploid_nucleotides[gt][0] == a && folded_diploid_nucleotides[gt][1] == b) ||
		folded_diploid_nucleotides[gt][0] == b && folded_diploid_nucleotides[gt][1] == a);
	int unfolded = a*4+b;
	int gt2 = folded_diploid_genotypes[unfolded];
	BOOST_CHECK((folded_diploid_nucleotides[gt2][0] == a && folded_diploid_nucleotides[gt2][1] == b) ||
		folded_diploid_nucleotides[gt2][0] == b && folded_diploid_nucleotides[gt2][1] == a);

	if(a <= b) {
		BOOST_CHECK(unfolded_diploid_genotypes_upper[gt2] == unfolded);
	}
	if(b <= a) {
		BOOST_CHECK(unfolded_diploid_genotypes_lower[gt2] == unfolded);
	}
}
