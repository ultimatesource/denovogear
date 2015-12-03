#define BOOST_TEST_MODULE "dng::utility::utility.h"

#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <utility>

#include "dng/utility.h"

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

  BOOST_CHECK(dng::utility::key_switch_tuple("", keys, keys[0]).second == 0);
  BOOST_CHECK(dng::utility::key_switch_tuple("0", keys, keys[0]).second == 0);
  BOOST_CHECK(dng::utility::key_switch_tuple("1", keys, keys[0]).second == 1);
  BOOST_CHECK(dng::utility::key_switch_tuple("2", keys, keys[0]).second == 2);
  BOOST_CHECK(dng::utility::key_switch_tuple("male", keys, keys[0]).second == 1);
  BOOST_CHECK(dng::utility::key_switch_tuple("female", keys, keys[0]).second == 2);
  BOOST_CHECK(dng::utility::key_switch_tuple("unknown", keys, keys[0]).second == 0);
  BOOST_CHECK(dng::utility::key_switch_tuple("MALE", keys, keys[0]).second == 1);
  BOOST_CHECK(dng::utility::key_switch_tuple("Female", keys, keys[0]).second == 2);

}


BOOST_AUTO_TEST_CASE(Test_parse_double_list)
{
  string nfreq = "0.3,0.2,0.2,0.3";
  vector<double> nfreq_a = {0.3, 0.2, 0.2, 0.3};
  char sep = ',';
  size_t size = 4;

  pair<vector<double>, bool> dup = dng::utility::parse_double_list(nfreq, sep, size);
  BOOST_CHECK(dup.second == true);
  BOOST_CHECK(dup.first.size() == size);
  BOOST_CHECK(dup.first == nfreq_a);

}

BOOST_AUTO_TEST_CASE(Test_to_pretty)
{
  BOOST_CHECK(dng::utility::to_pretty(1)=="1");
  BOOST_CHECK(dng::utility::to_pretty(1.234)=="1.234");
  BOOST_CHECK(dng::utility::to_pretty("string")=="string");
  BOOST_CHECK(dng::utility::to_pretty('c')=="c");
}

BOOST_AUTO_TEST_CASE(Test_extract_file_type)
{
  using namespace std;

  BOOST_CHECK(dng::utility::extract_file_type("foo.bar") == make_pair(string{"bar"},string{"foo.bar"}));
  BOOST_CHECK(dng::utility::extract_file_type("my:foo.bar") == make_pair(string{"my"},string{"foo.bar"}));
  BOOST_CHECK(dng::utility::extract_file_type(".bar") == make_pair(string{""},string{".bar"}));
  BOOST_CHECK(dng::utility::extract_file_type("my:.foo.bar") == make_pair(string{"my"},string{".foo.bar"}));
  BOOST_CHECK(dng::utility::extract_file_type(".foo.bar") == make_pair(string{"bar"},string{".foo.bar"}));
  BOOST_CHECK(dng::utility::extract_file_type("") == make_pair(string{""},string{""}));
  BOOST_CHECK(dng::utility::extract_file_type({}) == make_pair(string{},string{}));
  BOOST_CHECK(dng::utility::extract_file_type("C:\\foo.bar") == make_pair(string{"bar"},string{"C:\\foo.bar"}));
  BOOST_CHECK(dng::utility::extract_file_type("CC:C:\\foo.bar") == make_pair(string{"CC"},string{"C:\\foo.bar"}));

  BOOST_CHECK(dng::utility::extract_file_type(" \f\n\r\t\vfoo.bar \f\n\r\t\v") == make_pair(string{"bar"},string{"foo.bar"}));
  BOOST_CHECK(dng::utility::extract_file_type(" \f\n\r\t\vmy:foo.bar \f\n\r\t\v") == make_pair(string{"my"},string{"foo.bar"}));
  BOOST_CHECK(dng::utility::extract_file_type(" \f\n\r\t\v.bar \f\n\r\t\v") == make_pair(string{""},string{".bar"}));
  BOOST_CHECK(dng::utility::extract_file_type(" \f\n\r\t\v") == make_pair(string{},string{}));
}
