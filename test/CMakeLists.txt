#ADD_CUSTOM_TARGET(testdata
#  COMMAND "ENABLE_TESTING()"
#)


################################################################################
# Determine the location of the executables/scripts and any data sets required for full package testing
# TODO: Search for packages and data installed outside of denovogear directory, move to ../modules/
SET(DNG_DNM_EXE ${CMAKE_BINARY_DIR}/src/dng-dnm)
SET(DNG_PHASER_EXE ${CMAKE_BINARY_DIR}/src/dng-phaser)
SET(DNG_CALL_EXE ${CMAKE_BINARY_DIR}/src/dng-call)
IF(NOT TEST_DATA_DIR)
  IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/testdata)
    SET(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/testdata)
  ELSEIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/../testdata)
    SET(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/../testdata)
  ELSE(NOT TEST_DATA_DIR)
    MESSAGE(STATUS "WARNING: Unable to find data for full package testing. Please set TEST_DATA_DIR.")
  ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/testdata)
ENDIF(NOT TEST_DATA_DIR)


################################################################################
# Compile unit tests
include_directories(../src/include)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR}/unsupported)

IF(NOT Boost_USE_STATIC_LIBS)
  add_definitions(-DBOOST_TEST_DYN_LINK)
ENDIF()

# Make a test for each unit test file located in test/ dir
file(GLOB_RECURSE TESTS test_*[cc|cpp])
foreach(test ${TESTS})
  get_filename_component(testName ${test} NAME_WE)
  add_executable(${testName} EXCLUDE_FROM_ALL ${test})
  set_target_properties(${testName} PROPERTIES CXX_STANDARD 11)
  set_target_properties(${testName} PROPERTIES CXX_STANDARD_REQUIRED ON)
  target_link_libraries(${testName} Threads::Threads)
  IF(Boost_FOUND)
    target_link_libraries(${testName}
      ${Boost_PROGRAM_OPTIONS_LIBRARY}
      ${Boost_FILESYSTEM_LIBRARY}
      ${Boost_SYSTEM_LIBRARY}
      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    )
  ENDIF(Boost_FOUND)
  IF(HTSLIB_FOUND)
    target_link_libraries(${testName} ${HTSLIB_LIBRARIES})
  ENDIF(HTSLIB_FOUND)
  IF(ZLIB_FOUND)
    target_link_libraries(${testName} ZLIB::ZLIB)
  ENDIF(ZLIB_FOUND)
  add_test(${testName}_build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${testName})
  add_test(${testName}_run ${testName})
  set_tests_properties(${testName}_run PROPERTIES DEPENDS ${testName}_build)
endforeach(test)

################################################################################
# Full executable tests
# TODO: Currently can check dng found mutation, but not for false positives. Need stricter way to check output.
IF(DEFINED TEST_DATA_DIR AND NOT IS_DIRECTORY "${TEST_DATA_DIR}")
  MESSAGE(STATUS "WARNING: Could not find directory ${TEST_DATA_DIR}. Full Package tests will fail")
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

ADD_TEST(FullTest_SampleCEU_DNM ${DNG_DNM_EXE} auto --ped ${TEST_DATA_DIR}/sample_CEU/sample_CEU.ped --bcf ${TEST_DATA_DIR}/sample_CEU/sample_CEU.vcf --snp_mrate 2e-10 --indel_mrate 1e-11)
SET_TESTS_PROPERTIES(FullTest_SampleCEU_DNM PROPERTIES
  PASS_REGULAR_EXPRESSION "DENOVO-SNP CHILD_ID: NA12878_vald-sorted\.bam\.bam" 
  PASS_REGULAR_EXPRESSION "chr: 2"
  PASS_REGULAR_EXPRESSION "pos: 214668360" 
  PASS_REGULAR_EXPRESSION "ref: G"
  PASS_REGULAR_EXPRESSION "alt: A"
  PASS_REGULAR_EXPRESSION "maxlike_null: 3\.95324e-12" 
  PASS_REGULAR_EXPRESSION "pp_null: 0\.0195898"
  PASS_REGULAR_EXPRESSION "tgt_null\(child\/mom\/dad\): GG\/GG\/GG" 
  PASS_REGULAR_EXPRESSION "snpcode: 1" 
  PASS_REGULAR_EXPRESSION "code: 6"
  PASS_REGULAR_EXPRESSION "maxlike_dnm: 1\.98602e-10" 
  PASS_REGULAR_EXPRESSION "pp_dnm: 0\.98041"
  PASS_REGULAR_EXPRESSION "tgt_dnm\(child\/mom\/dad\): AG\/GG\/GG" 
  PASS_REGULAR_EXPRESSION "lookup: 4" 
  PASS_REGULAR_EXPRESSION "flag: 0"
  PASS_REGULAR_EXPRESSION "READ_DEPTH child: 48 dad: 76 mom: 34"
  PASS_REGULAR_EXPRESSION "MAPPING_QUALITY child: 59 dad: 59 mom: 59")

ADD_TEST(FullTest_SamplePaired_DNM ${DNG_DNM_EXE} auto --ped ${TEST_DATA_DIR}/sample_CEU/sample_paired.ped --bcf ${TEST_DATA_DIR}/sample_CEU/sample_CEU.vcf)
SET_TESTS_PROPERTIES(FullTest_SamplePaired_DNM PROPERTIES
  PASS_REGULAR_EXPRESSION "TUMOR_ID: NA12878_vald-sorted\.bam\.bam"
  PASS_REGULAR_EXPRESSION "NORMAL_ID: NA12891_vald-sorted\.bam\.bam"
  PASS_REGULAR_EXPRESSION "chr: 2"
  PASS_REGULAR_EXPRESSION "pos: 214668360" 
  PASS_REGULAR_EXPRESSION "ref: G" 
  PASS_REGULAR_EXPRESSION "alt: A" 
  PASS_REGULAR_EXPRESSION "maxlike_null: 3\.98107e-12" 
  PASS_REGULAR_EXPRESSION "pp_null: 0\.00396529"
  PASS_REGULAR_EXPRESSION "tgt_null\(normal\/tumor\): GG\/GG" 
  PASS_REGULAR_EXPRESSION "maxlike_dnm: 1e-09" 
  PASS_REGULAR_EXPRESSION "pp_dnm: 0\.996035"
  PASS_REGULAR_EXPRESSION "tgt_dnm\(normal\/tumor\): GG\/AG"
  PASS_REGULAR_EXPRESSION "READ_DEPTH tumor: 48"
  PASS_REGULAR_EXPRESSION "normal: 76"
  PASS_REGULAR_EXPRESSION "tumor: 59" 
  PASS_REGULAR_EXPRESSION "normal: 59"
  PASS_REGULAR_EXPRESSION "null_snpcode: 1" 
  PASS_REGULAR_EXPRESSION "dnm_snpcode: 2"
)

ADD_TEST(FullTest_Sample_Phaser ${DNG_PHASER_EXE} --dnm ${TEST_DATA_DIR}/sample_Phaser/sample_phasing_dnm_f --pgt ${TEST_DATA_DIR}/sample_Phaser/sample_phasing_GTs_f --bam ${TEST_DATA_DIR}/sample_Phaser/test1.bam)
SET_TESTS_PROPERTIES(FullTest_Sample_Phaser PROPERTIES
  PASS_REGULAR_EXPRESSION "DNM_pos 1:75884343	INHERITED T	VARIANT C"
  PASS_REGULAR_EXPRESSION "HAP POS 75884200 p1: AT p2: AC"
  PASS_REGULAR_EXPRESSION "Base at DNM position: A Base at phasing position: A	 INFERRED PARENT OF ORIGIN for DNM: N/A SUPPORTING READ COUNT: 1"
  PASS_REGULAR_EXPRESSION "Base at DNM position: G Base at phasing position: C	 INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 1"
  PASS_REGULAR_EXPRESSION "DNM_pos 1:110583335	INHERITED G	VARIANT A"
  PASS_REGULAR_EXPRESSION "DNM_pos 1:182974758	INHERITED G	VARIANT A"
  PASS_REGULAR_EXPRESSION "HAP POS 182974760 p1: GT p2: CA"
  PASS_REGULAR_EXPRESSION "Base at DNM position: G Base at phasing position: G	 INFERRED PARENT OF ORIGIN for DNM: p2 SUPPORTING READ COUNT: 1"
)


