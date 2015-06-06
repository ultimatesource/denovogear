#define BOOST_TEST_MODULE "hts::extra.h"

#include <boost/test/unit_test.hpp>

#include "dng/hts/extra.h"

BOOST_AUTO_TEST_CASE(Test_extract_file_type)
{
  std::pair<std::string, std::string> duple;

  duple = hts::extra::extract_file_type("bcf:test1.bcf");
  BOOST_CHECK(duple.first.compare("bcf")==0);
  BOOST_CHECK(duple.second.compare("test1.bcf")==0);

  duple = hts::extra::extract_file_type("vcf:test2.vcf");
  BOOST_CHECK(duple.first.compare("vcf")==0);
  BOOST_CHECK(duple.second.compare("test2.vcf")==0);

  duple = hts::extra::extract_file_type("bcf:test3");
  BOOST_CHECK(duple.first.compare("bcf")==0);
  BOOST_CHECK(duple.second.compare("test3")==0);

  duple = hts::extra::extract_file_type("vcf:test4");
  BOOST_CHECK(duple.first.compare("vcf")==0);
  BOOST_CHECK(duple.second.compare("test4")==0);

  duple = hts::extra::extract_file_type("test5.bcf");
  BOOST_CHECK(duple.first.compare("bcf")==0);
  BOOST_CHECK(duple.second.compare("test5.bcf")==0);

  duple = hts::extra::extract_file_type("test6.vcf");
  BOOST_CHECK(duple.first.compare("vcf")==0);
  BOOST_CHECK(duple.second.compare("test6.vcf")==0);

  duple = hts::extra::extract_file_type("test7");
  BOOST_CHECK(duple.first.compare("")==0);
  BOOST_CHECK(duple.second.compare("test7")==0);

}
