#define BOOST_TEST_MODULE dng::ReadGroups

#include <vector>
#include <string>
#include <utility>

#include <dng/read_group.h>

using namespace std;


BOOST_AUTO_TEST_CASE(test_parse_header_text) {
    const char buffer[] =
        "@RG\tID:a\tLB:a\tSM:a\n"
        "@RG\tID:b\tLB:b\tSM:b\n"
        "@RG\tID:c1\tLB:c\tSM:c\n"
        "@RG\tID:c2\tLB:c\tSM:c\n"
        ;
    dng::ReadGroups rgs;
    rgs.ParseHeaderText(buffer, "LB");
    BOOST_CHECK(rgs.library_from_index(0) == 0);
    BOOST_CHECK(rgs.library_from_index(1) == 1);
    BOOST_CHECK(rgs.library_from_index(2) == 2);
    BOOST_CHECK(rgs.library_from_index(3) == 2);

    vector<string> bad_libs = {"b"};
    rgs.EraseLibraries(bad_libs);

    BOOST_CHECK(rgs.library_from_index(0) == 0);
    BOOST_CHECK(rgs.library_from_index(1) == 1);
    BOOST_CHECK(rgs.library_from_index(2) == 1);
}


