#define BOOST_TEST_MODULE "dng::util::seq.h"

#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>

#include "dng/seq.h"

using namespace std;

BOOST_AUTO_TEST_CASE(Test_char_index)
{
  BOOST_CHECK(dng::seq::char_index('N') == 4);
  BOOST_CHECK(dng::seq::char_index('A') == 0);
  BOOST_CHECK(dng::seq::char_index('C') == 1);
  BOOST_CHECK(dng::seq::char_index('G') == 2);
  BOOST_CHECK(dng::seq::char_index('T') == 3);
}

BOOST_AUTO_TEST_CASE(Test_indexed_char)
{
  BOOST_CHECK(dng::seq::indexed_char(4) == 'N');
  BOOST_CHECK(dng::seq::indexed_char(0) == 'A');
  BOOST_CHECK(dng::seq::indexed_char(1) == 'C');
  BOOST_CHECK(dng::seq::indexed_char(2) == 'G');
  BOOST_CHECK(dng::seq::indexed_char(3) == 'T');

}
