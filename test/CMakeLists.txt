FIND_PACKAGE(Boost COMPONENTS unit_test_framework REQUIRED)
FIND_PACKAGE(HTSLIB)

################################################################################
# Determine the location of the executables/scripts and any data sets required for full package testing
# TODO: Search for packages and data installed outside of denovogear directory, move to ../modules/
SET(DNG_CALL_EXE ${CMAKE_BINARY_DIR}/src/dng-call)
IF(NOT TEST_DATA_DIR)
  IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/testdata)
    SET(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/testdata)
  ELSEIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/../testdata)
    SET(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/../testdata)
  ELSE(NOT TEST_DATA_DIR)
    MESSAGE("Unable to find data for full package testing. Please set TEST_DATA_DIR")
  ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/testdata)
ENDIF(NOT TEST_DATA_DIR)
SET(TEST_DATA_DIR /home/kdai1/workspace/testdata)


################################################################################
# Compile unit tests
include_directories(../src/include)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR}/unsupported)
add_definitions(-DBOOST_TEST_DYN_LINK)

# Make a test for each unit test file located in test/ dir
file(GLOB_RECURSE TESTS test_*[cc|cpp])
foreach(test ${TESTS})
  get_filename_component(testName ${test} NAME_WE)
  add_executable(${testName} ${test})
  set_target_properties(${testName} PROPERTIES CXX_STANDARD 11)
  set_target_properties(${testName} PROPERTIES CXX_STANDARD_REQUIRED ON)
  target_link_libraries(${testName} ${Boost_LIBRARIES} ${HTSLIB_LIBRARIES})
  add_test(${testName} ${testName})
endforeach(test)

################################################################################
# Full executable tests
IF(NOT IS_DIRECTORY "${TEST_DATA_DIR}")
  MESSAGE("Could not find directory ${TEST_DATA_DIR}. Full Package tests will fail")
ENDIF()
    
ADD_TEST(Test_Help ${DNG_CALL_EXE})
SET_TESTS_PROPERTIES(Test_Help PROPERTIES PASS_REGULAR_EXPRESSION "Usage:")

ADD_TEST(FullTest_Sample5.3_BAM ${DNG_CALL_EXE} -f ${TEST_DATA_DIR}/sample_5_3/sample-5.3_ref.fasta.gz -p ${TEST_DATA_DIR}/sample_5_3/ceu.ped -m 0.001 ${TEST_DATA_DIR}/sample_5_3/test1.bam)
SET_TESTS_PROPERTIES(FullTest_Sample5.3_BAM PROPERTIES PASS_REGULAR_EXPRESSION "^.*(126385924)")

ADD_TEST(FullTest_Sample5.3_VCF ${DNG_CALL_EXE} -p ${TEST_DATA_DIR}/sample_5_3/ceu.ped -m 0.001 ${TEST_DATA_DIR}/sample_5_3/test1.vcf)
SET_TESTS_PROPERTIES(FullTest_Sample5.3_VCF PROPERTIES 
  PASS_REGULAR_EXPRESSION "LL=-74\.68" 
  PASS_REGULAR_EXPRESSION "25,29,1"
  PASS_REGULAR_EXPRESSION "57,0,0"
  PASS_REGULAR_EXPRESSION "76,1,0")
