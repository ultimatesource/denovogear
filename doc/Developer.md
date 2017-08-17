# Developer Information

## Compile

```bash
mkdir ./build/clangdebug
CXX=clang++ CC=clang cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug
ninja
ninja test
```

## Measure Code Coverage

```bash
mkdir ./build/coverage
CXX=g++ CC=gcc cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DDNG_DEVEL_ENABLE_COVERAGE_REPORT=1
ninja
```
### Run Unit Tests only

```bash
ctest --output-on-failure -R Unit
gcovr -d -r . --html --html-details -o output.html
```

### Remove Coverage Files

```bash
find . -name '*.gcda' -delete
find . -name '*.gcov' -delete
```

# Writing Unit Tests

Unit tests for the `dng` namespace are placed in the `Tests/Unit/dng` folder and registered
by editing the end of `Tests/Unit/CMakeLists.txt`. Unit testing is done via the [Boost
unit-testing framework](http://www.boost.org/doc/libs/1_64_0/libs/test/doc/html/index.html).

## Example

```c++
#define BOOST_TEST_MODULE dng::utility
#include <dng/utility.h>
#include "../testing.h"

using namespace dng::utility;

BOOST_AUTO_TEST_CASE(test_locations) {
    BOOST_CHECK_EQUAL(LOCATION_MAX, 0x7FFFFFFF7FFFFFFFll);

    auto test = [](int t, int p) -> void {
    BOOST_TEST_CONTEXT("t=" << t << ", p=" << p) {
        BOOST_REQUIRE_GE(t, 0);
        BOOST_REQUIRE_GE(p, 0);
        BOOST_REQUIRE_LE(t, INT32_MAX);
        BOOST_REQUIRE_LE(p, INT32_MAX);

        location_t test_location = make_location(t,p);
        int64_t expected_location = ((int64_t)t << 32) | p;
        BOOST_CHECK_EQUAL(test_location, expected_location);
        BOOST_CHECK_EQUAL(location_to_contig(test_location), t);
        BOOST_CHECK_EQUAL(location_to_position(test_location), p);
    }};

    test(0,0);
    test(0,1);
    test(10,2147483647);
    test(2147483647,2147483647);
}
```

## Tips

 - A unit test file should correspond to either a namespace, class, or source-code file.
 - A unit test must be self-contained. If a test needs to access data, embed that data in the
   source file and either read it from memory or write it to disk using `dng::detail::AutoTempFile`.
   If the data is too large for this, use an integration test and the testdata repo.
 - Put the header of the namespace or class being tested first. This will catch any missing includes.
 - Use lambda functions to keep testing logic together with test-cases.
 - Use `Test/Unit/xorshift64.h` to generate random testing data.
 - Unit tests must support both Boost.Test 3 (primarily) and 2 (secondarily).
   (Some distros still ship with Boost.Test 2.)
   `Tests/Unit/testing.h` contains compatibility macros. Create new ones as needed.

## Accessing Private Members

Some unit tests are made easier by accessing private members of a class.
Use the `DNG_UNIT_TEST_CLASS` and `DNG_UNIT_TEST_FUNCTION` to declare friend
classes and functions that can access private and protected members.

```c++
// src/include/dng/example.h

#include <dng/detail/unit_test.h>

namespace dng {
class Example {
    Example(double d) : d_{d} { }
private:
    double d_;
    DNG_UNIT_TEST_CLASS(unittest_dng_example);
};

}
```

```c++
// Tests/Unit/dng/example.cc

#include <dng/example.h>
#include "../testing.h"

namespace dng {
    struct unittest_dng_example {
    GETTERS_FOR_MEMBER_VARIABLES(RelationshipGraph,
        (d)
    };
}

using u = dng::unittest_dng_example;
using namespace dng;

BOOST_AUTO_TEST_CASE(test_Example_construct) {
    Example test(10.0);
    BOOST_CHECK_EQUAL(u::d(test), 10.0);
}

```
