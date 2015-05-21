#define BOOST_TEST_MODULE "dng::util::utilities.h"

#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>

#include "dng/utilities.h"

using namespace std;


BOOST_AUTO_TEST_CASE(Test_key_switch_tuple)
{
  pair<string,int> keys[] = {
    {"0", 0},
    {"1", 1},
    {"2", 2},
    {"male", 1},
    {"female", 2},
    {"unknown", 0}
  };

  BOOST_CHECK(dng::util::key_switch_tuple("", keys, keys[0]).second == 0);
  BOOST_CHECK(dng::util::key_switch_tuple("0", keys, keys[0]).second == 0);
  BOOST_CHECK(dng::util::key_switch_tuple("1", keys, keys[0]).second == 1);
  BOOST_CHECK(dng::util::key_switch_tuple("2", keys, keys[0]).second == 2);
  BOOST_CHECK(dng::util::key_switch_tuple("male", keys, keys[0]).second == 1);
  BOOST_CHECK(dng::util::key_switch_tuple("female", keys, keys[0]).second == 2);
  BOOST_CHECK(dng::util::key_switch_tuple("unknown", keys, keys[0]).second == 0);
  BOOST_CHECK(dng::util::key_switch_tuple("MALE", keys, keys[0]).second == 1);
  BOOST_CHECK(dng::util::key_switch_tuple("Female", keys, keys[0]).second == 2);

}


BOOST_AUTO_TEST_CASE(Test_parse_double_list)
{
  string nfreq = "0.3,0.2,0.2,0.3";
  vector<double> nfreq_a = {0.3, 0.2, 0.2, 0.3};
  char sep = ',';
  size_t size = 4;

  pair<vector<double>, bool> dup = dng::util::parse_double_list(nfreq, sep, size);
  BOOST_CHECK(dup.second == true);
  BOOST_CHECK(dup.first.size() == size);
  BOOST_CHECK(dup.first == nfreq_a);

}

BOOST_AUTO_TEST_CASE(Test_to_pretty)
{
  BOOST_CHECK(dng::util::to_pretty(1)=="1");
  BOOST_CHECK(dng::util::to_pretty(1.234)=="1.234");
  BOOST_CHECK(dng::util::to_pretty("string")=="string");
  BOOST_CHECK(dng::util::to_pretty('c')=="c");
}
