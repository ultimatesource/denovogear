
################################################################################
# Determine the location of the executables/scripts and any data sets required for full package testing
# TODO: Search for packages and data installed outside of denovogear directory, move to ../modules/
SET(DNG_DNM_EXE ${CMAKE_BINARY_DIR}/src/dng-dnm)
SET(DNG_PHASER_EXE ${CMAKE_BINARY_DIR}/src/dng-phaser)
SET(DNG_CALL_EXE ${CMAKE_BINARY_DIR}/src/dng-call)
SET(DNG_DNM_EXE ${CMAKE_BINARY_DIR}/src/dng-dnm)
SET(DNG_PHASER_EXE ${CMAKE_BINARY_DIR}/src/dng-phaser)
IF(NOT TEST_DATA_DIR)
  IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/testdata)
    SET(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/testdata)
  ELSEIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/../testdata)
    SET(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/../testdata)
#  ELSE(NOT TEST_DATA_DIR)
#    MESSAGE(STATUS "WARNING: Unable to find data for full package testing. Please set TEST_DATA_DIR.")
  ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/testdata)
ENDIF(NOT TEST_DATA_DIR)


###############################################################################
# Tests on the source code
ADD_TEST(build_gitkeep ${CMAKE_COMMAND} -E md5sum "${CMAKE_SOURCE_DIR}/build/.gitkeep")
SET_TESTS_PROPERTIES(build_gitkeep PROPERTIES
  PASS_REGULAR_EXPRESSION "d41d8cd98f00b204e9800998ecf8427e"
)



################################################################################
# Compile unit tests
include_directories(../src/include)


# Make a test for each unit test file located in test/ dir
file(GLOB_RECURSE TESTS test_*[cc|cpp])
foreach(test ${TESTS})
  get_filename_component(testName ${test} NAME_WE)
  add_executable(${testName} EXCLUDE_FROM_ALL ${test})
  target_link_libraries(${testName}
    Threads::Threads
    Boost::PROGRAM_OPTIONS
    Boost::FILESYSTEM
    Boost::SYSTEM
    Boost::UNIT_TEST_FRAMEWORK
    EIGEN3::EIGEN3
    HTSLIB::HTSLIB
  )
  add_test(${testName}_build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${testName})
  add_test(${testName}_run ${testName})
  set_tests_properties(${testName}_run PROPERTIES DEPENDS ${testName}_build)
endforeach(test)

################################################################################



################################################################################
## TODO: How to build something like "make timetrail" and compile all time measuring files
IF(DEVEL_MODE)
  ADD_DEFINITIONS(-DDNG_DEVEL)
ENDIF()
file(GLOB_RECURSE TESTS2 "testt_time*[cc|cpp]")
foreach(test ${TESTS2})
  get_filename_component(testName ${test} NAME_WE)
  add_executable(${testName} EXCLUDE_FROM_ALL ${test}
          ${CMAKE_SOURCE_DIR}/src/lib/peeling.cc
          ##            ../src/lib/likelihood.cc ../src/lib/newick.cc ../src/lib/pedigree.cc
          ##            ../src/lib/mutation.cc ../src/lib/stats.cc
          )
  target_link_libraries(${testName}
          Threads::Threads
          Boost::PROGRAM_OPTIONS
          Boost::FILESYSTEM
          Boost::SYSTEM
          Boost::UNIT_TEST_FRAMEWORK
          EIGEN3::EIGEN3
          HTSLIB::HTSLIB
#          Boost::TIMER
          )
#  add_test(${testName}_build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${testName})
#  add_test(${testName}_run ${testName})
#  set_tests_properties(${testName}_run PROPERTIES DEPENDS ${testName}_build)
endforeach(test)





## temp Compile unit tests with *.cc
# TODO: Use a different filename pattern for now (test2_* instead of test_*).
# TODO: How to add min *.cc files instead of all *.cc later?
include_directories(./)
include_directories(./lib/)

file(GLOB_RECURSE TESTS2 "test2_*[cc|cpp]")
foreach(test ${TESTS2})
  get_filename_component(testName ${test} NAME_WE)
  add_executable(${testName} EXCLUDE_FROM_ALL ${test}
          ${CMAKE_SOURCE_DIR}/src/lib/peeling.cc
          ##            ../src/lib/likelihood.cc ../src/lib/newick.cc ../src/lib/pedigree.cc
          ##            ../src/lib/mutation.cc ../src/lib/stats.cc
          )
  target_link_libraries(${testName}
          Threads::Threads
          Boost::PROGRAM_OPTIONS
          Boost::FILESYSTEM
          Boost::SYSTEM
          Boost::UNIT_TEST_FRAMEWORK
          EIGEN3::EIGEN3
          HTSLIB::HTSLIB
          )
  add_test(${testName}_build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${testName})
  add_test(${testName}_run ${testName})
  set_tests_properties(${testName}_run PROPERTIES DEPENDS ${testName}_build)
endforeach(test)

################################################################################


# Full executable tests

# Check if using download source package or running latest code inside a git repository
IF(NOT DEFINED DNG_VERSION_COMMIT OR DNG_VERSION_COMMIT STREQUAL "unknown")
  SET(DATA_PACKAGE "v${DNG_VERSION_MAJOR}.${DNG_VERSION_MINOR}")
  IF(DEFINED DNG_VERSION_PATCH)
    SET(DATA_PACKAGE "${DATA_PACKAGE}.${DNG_VERSION_PATCH}")
  ENDIF()
  SET(DATA_PACKAGE "${DATA_PACKAGE}.tar.gz")
  SET(DOWNLOAD_URL "https://github.com/kaeldai/testdata/archive/${DATA_PACKAGE}")
ENDIF()

IF(NOT DEFINED TEST_DATA_DIR)
  # set test data download directory if not already set
  SET(TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/testdata")
ENDIF()

IF(DEFINED DOWNLOAD_URL)
  # if user testing code from package then get the corresponding archived tests.
  ExternalProject_Add(testdata
    DEPENDS dng-phaser dng-dnm dng-call
    PREFIX "testdata"
    EXCLUDE_FROM_ALL 1
    DOWNLOAD_DIR "${TEST_DATA_DIR}"
    URL "https://github.com/kaeldai/testdata/archive/${DATA_PACKAGE}"
    SOURCE_DIR "${TEST_DATA_DIR}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")
ELSE()
  # if code was downloaded directly from github, then grab the lastest tests
  ExternalProject_Add(testdata
    DEPENDS dng-phaser dng-dnm dng-call
    PREFIX "testdata"
    EXCLUDE_FROM_ALL TRUE
    DOWNLOAD_DIR "${TEST_DATA_DIR}"
    GIT_REPOSITORY "https://github.com/kaeldai/testdata.git"
    SOURCE_DIR "${TEST_DATA_DIR}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")
ENDIF()

IF(NOT DEFINED TEST_DATA_DIR OR (DEFINED TEST_DATA_DIR AND NOT IS_DIRECTORY "${TEST_DATA_DIR}"))
  MESSAGE(STATUS "WARNING: Could not find directory ${TEST_DATA_DIR}. Full package tests will fail.")
  MESSAGE(STATUS "         Please run \'cmake testdata\' before running \'make test\'")
ENDIF()
    
ADD_TEST(CALL_Help ${DNG_CALL_EXE})
SET_TESTS_PROPERTIES(CALL_Help PROPERTIES PASS_REGULAR_EXPRESSION "Usage:")

ADD_TEST(CALL_Sample5.3_BAM ${DNG_CALL_EXE} -f ${TEST_DATA_DIR}/sample_5_3/sample-5.3_ref.fasta.gz -p ${TEST_DATA_DIR}/sample_5_3/ceu.ped -m 0.001 ${TEST_DATA_DIR}/sample_5_3/test1.bam)
SET_TESTS_PROPERTIES(CALL_Sample5.3_BAM PROPERTIES PASS_REGULAR_EXPRESSION "^.*(126385924)")
SET_TESTS_PROPERTIES(CALL_Sample5.3_BAM PROPERTIES 
  PASS_REGULAR_EXPRESSION "^.*(126385924)"
)

ADD_TEST(CALL_Sample5.3_VCF ${DNG_CALL_EXE} -p ${TEST_DATA_DIR}/sample_5_3/ceu.ped -m 0.001 ${TEST_DATA_DIR}/sample_5_3/test1.vcf)
SET_TESTS_PROPERTIES(CALL_Sample5.3_VCF PROPERTIES 
  PASS_REGULAR_EXPRESSION "LL=-74\.68" 
  PASS_REGULAR_EXPRESSION "25,29,1"
  PASS_REGULAR_EXPRESSION "57,0,0"
  PASS_REGULAR_EXPRESSION "76,1,0"
  PASS_REGULAR_EXPRESSION "MUP=0\.118233"
  PASS_REGULAR_EXPRESSION "LLD=-30\.5955"
  PASS_REGULAR_EXPRESSION "LLH=-6\.68771" 
  PASS_REGULAR_EXPRESSION "MUX=0\.118233" 
  PASS_REGULAR_EXPRESSION "MU1P=0\.118233"
  PASS_REGULAR_EXPRESSION "DNT=GGxGG>GT"
  PASS_REGULAR_EXPRESSION "DNL=LB-NA12878:Solexa-135852"
  PASS_REGULAR_EXPRESSION "DNQ=42"
  PASS_REGULAR_EXPRESSION "DNC=100"
  PASS_REGULAR_EXPRESSION "DP=189"
  PASS_REGULAR_EXPRESSION "AD=158,30,1"
  PASS_REGULAR_EXPRESSION "0/0:38:0\.999831,0\.000169293,1\.81217e-12,2\.76783e-11,3\.62326e-16,4\.77512e-20"
  PASS_REGULAR_EXPRESSION "0/0:18:0\.984381,0\.0156192,1\.68954e-10,1\.05777e-11,1\.69591e-14,1\.81182e-20"  
  PASS_REGULAR_EXPRESSION "0/0:9:0\.865987,0\.134013,1\.23154e-15,2\.16126e-09,4\.20711e-17,6\.32158e-26:-6\.74916,0,-5\.86931,-7\.56277,-7\.41299,-15\.7054:55:25,29,1:0\.999933:0\.999933"
  PASS_REGULAR_EXPRESSION "0/0:38:0\.999835,0\.000165352,3\.10306e-15,2\.79305e-11,6\.20477e-19,8\.17697e-23:0,-6\.64431,-17\.371,-6\.64431,-17\.371,-17\.371:57:57,0,0:3\.33278e-05:3\.33278e-05"
  PASS_REGULAR_EXPRESSION "0/0:18:0\.984385,0\.0156153,1\.0244e-13,1\.02758e-11,1\.0262e-17,5\.33554e-24:0,-4\.66917,-15\.8523,-7\.07309,-16\.1525,-18\.5564:77:76,1,0:3\.33292e-05:3\.33292e-05"
)

ADD_TEST(DNM_Help ${DNG_DNM_EXE} --h)
SET_TESTS_PROPERTIES(DNM_Help PROPERTIES PASS_REGULAR_EXPRESSION "Usage:")

ADD_TEST(DNM_InputErr ${DNG_DNM_EXE} auto)
SET_TESTS_PROPERTIES(DNM_InputErr PROPERTIES PASS_REGULAR_EXPRESSION "ERROR")

ADD_TEST(DNM_SampleCEU ${DNG_DNM_EXE} auto --ped ${TEST_DATA_DIR}/sample_CEU/sample_CEU.ped --bcf ${TEST_DATA_DIR}/sample_CEU/sample_CEU.vcf --snp_mrate 2e-10 --indel_mrate 1e-11)
SET_TESTS_PROPERTIES(DNM_SampleCEU PROPERTIES
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
  PASS_REGULAR_EXPRESSION "MAPPING_QUALITY child: 59 dad: 59 mom: 59"
)

ADD_TEST(DNM_SamplePaired ${DNG_DNM_EXE} auto --ped ${TEST_DATA_DIR}/sample_CEU/sample_paired.ped --bcf ${TEST_DATA_DIR}/sample_CEU/sample_CEU.vcf)
SET_TESTS_PROPERTIES(DNM_SamplePaired PROPERTIES
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

ADD_TEST(DNM_ATMutation ${DNG_DNM_EXE} auto --ped ${TEST_DATA_DIR}/trio/mutationAT.ped --vcf ${TEST_DATA_DIR}/trio/mutationAT.vcf --snp_mrate 1e-8 --poly_rate 0.001)
SET_TESTS_PROPERTIES(DNM_ATMutation PROPERTIES
  PASS_REGULAR_EXPRESSION "chr: 2"
  PASS_REGULAR_EXPRESSION "pos: 214668361"
  PASS_REGULAR_EXPRESSION "maxlike_null: 9\.9301e-11" 
  PASS_REGULAR_EXPRESSION "pp_null: 0\.00990728"
  PASS_REGULAR_EXPRESSION "tgt_null(child\/mom\/dad): AA\/AA\/AA"
  PASS_REGULAR_EXPRESSION "snpcode: 1"
  PASS_REGULAR_EXPRESSION "code: 6"
  PASS_REGULAR_EXPRESSION "maxlike_dnm: 9\.9301e-09" 
  PASS_REGULAR_EXPRESSION "pp_dnm: 0\.990093"
  PASS_REGULAR_EXPRESSION "tgt_dnm(child\/mom\/dad): AT\/AA\/AA"
  PASS_REGULAR_EXPRESSION "lookup: 5" 
  PASS_REGULAR_EXPRESSION "flag: 0" 
)


# Test the help command 
ADD_TEST(Phaser_Help ${DNG_PHASER_EXE} --h)
SET_TESTS_PROPERTIES(Phaser_Help PROPERTIES PASS_REGULAR_EXPRESSION "Usage:")

# Test for input error message
ADD_TEST(Phaser_InputErr ${DNG_PHASER_EXE})
SET_TESTS_PROPERTIES(Phaser_InputErr PROPERTIES PASS_REGULAR_EXPRESSION "INPUT ERROR!")

ADD_TEST(Phaser_Sample ${DNG_PHASER_EXE} --dnm ${TEST_DATA_DIR}/sample_Phaser/sample_phasing_dnm_f --pgt ${TEST_DATA_DIR}/sample_Phaser/sample_phasing_GTs_f --bam ${TEST_DATA_DIR}/sample_Phaser/test1.bam)
SET_TESTS_PROPERTIES(Phaser_Sample PROPERTIES
  PASS_REGULAR_EXPRESSION "DNM_pos 1:75884343	INHERITED T	VARIANT C"
  PASS_REGULAR_EXPRESSION "HAP POS 75884200 p1: AT p2: AC"
  PASS_REGULAR_EXPRESSION "Base at DNM position: A Base at phasing position: A	 INFERRED PARENT OF ORIGIN for DNM: N/A SUPPORTING READ COUNT: 1"
  PASS_REGULAR_EXPRESSION "Base at DNM position: G Base at phasing position: C	 INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 1"
  PASS_REGULAR_EXPRESSION "DNM_pos 1:110583335	INHERITED G	VARIANT A"
  PASS_REGULAR_EXPRESSION "DNM_pos 1:182974758	INHERITED G	VARIANT A"
  PASS_REGULAR_EXPRESSION "HAP POS 182974760 p1: GT p2: CA"
  PASS_REGULAR_EXPRESSION "Base at DNM position: G Base at phasing position: G	 INFERRED PARENT OF ORIGIN for DNM: p2 SUPPORTING READ COUNT: 1"
)

##### Issue 51 (and 52) - Test cases for various library identifiers (by LB, SM, and RG). ####
# Only check the headers. Data will be different (and possibly meaningless) depending on how read-group libraries are grouped together 
# TODO: Incorporate test into the test_read_groups.cc unit-tests.
ADD_TEST(CALL_Sample5.3_BAM_LB_TAG ${DNG_CALL_EXE} -f ${TEST_DATA_DIR}/sample_5_3/sample-5.3_ref.fasta.gz -p ${TEST_DATA_DIR}/sample_5_3/ceu.ped -m 0.001 --rgtag "LB" ${TEST_DATA_DIR}/sample_5_3/test1.bam)
SET_TESTS_PROPERTIES(CALL_Sample5.3_BAM_LB_TAG PROPERTIES 
	PASS_REGULAR_EXPRESSION "FORMAT\tGL-1\tGL-2\tLB-NA12878-Solexa-135852\tLB-NA12891-Solexa-135851\tLB-NA12892-Solexa-135853"
)
ADD_TEST(CALL_Sample5.3_BAM_SM_TAG ${DNG_CALL_EXE} -f ${TEST_DATA_DIR}/sample_5_3/sample-5.3_ref.fasta.gz -p ${TEST_DATA_DIR}/sample_5_3/ceu.ped -m 0.001 --rgtag "SM" ${TEST_DATA_DIR}/sample_5_3/test1.bam)
SET_TESTS_PROPERTIES(CALL_Sample5.3_BAM_SM_TAG PROPERTIES 
	PASS_REGULAR_EXPRESSION "FORMAT\tGL-1\tGL-2\tLB-NA12878\tLB-NA12891\tLB-NA12892"
)
ADD_TEST(CALL_Sample5.3_BAM_ID_TAG ${DNG_CALL_EXE} -f ${TEST_DATA_DIR}/sample_5_3/sample-5.3_ref.fasta.gz -p ${TEST_DATA_DIR}/sample_5_3/ceu.ped -m 0.001 --rgtag "ID" ${TEST_DATA_DIR}/sample_5_3/test1.bam)
SET_TESTS_PROPERTIES(CALL_Sample5.3_BAM_ID_TAG PROPERTIES 
	PASS_REGULAR_EXPRESSION "FORMAT\tGL-1\tGL-2\tSM-NA12891\tSM-NA12892\tSM-NA12878\tLB-NA12878-H06HD.1\tLB-NA12878-H06HD.2\tLB-NA12878-H06JU.1\tLB-NA12891-H03N7.1\tLB-NA12891-H03N7.2\tLB-NA12891-H05F1.2\tLB-NA12892-H06JH.1\tLB-NA12892-H06JH.2\tLB-NA12892-H06JU.2"
)

##### Issue 55 - Separate header and bam file
ADD_TEST(CALL_Sample5.3_SAM_SEP_HDR ${DNG_CALL_EXE} -f ${TEST_DATA_DIR}/sep_header/sample-5.3_ref.fasta.gz -p ${TEST_DATA_DIR}/sep_header/ceu.ped -m 0.001 -h ${TEST_DATA_DIR}/sep_header/test1_hdr.sam ${TEST_DATA_DIR}/sep_header/test1_nohdr.sam)
SET_TESTS_PROPERTIES(CALL_Sample5.3_SAM_SEP_HDR PROPERTIES 
	PASS_REGULAR_EXPRESSION "FORMAT\tGL-1\tGL-2\tLB-NA12878-Solexa-135852\tLB-NA12891-Solexa-135851\tLB-NA12892-Solexa-135853"
)
