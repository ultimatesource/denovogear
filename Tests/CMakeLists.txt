###############################################################################
# Testdata information
set(TESTDATA_COMMIT b4e91b8042a2b265131d9985d33f43e55f457a94)
set(TESTDATA_SHA256 470396a0924823dfa76fc53785348dc56e249a9889fcbaddacf58b7912857e6a)
set(TESTDATA_REPO denovogear/testdata)
set(TESTDATA_URL "https://api.github.com/repos/${TESTDATA_REPO}/tarball/${TESTDATA_COMMIT}")
set(TESTDATA_DIR "${CMAKE_BINARY_DIR}/testdata")

# Set TESTDATA to specify a directory (e.g. git clone) that contains the testing data
set(TESTDATA "" CACHE PATH "Source directory of custom testing data.")


# If TESTDATA is not set download the testing data
if(NOT TESTDATA)
  set(testdata_args
    DOWNLOAD_NAME testdata.tar.gz
    DOWNLOAD_DIR "${CMAKE_BINARY_DIR}"
    URL "${TESTDATA_URL}"
    URL_HASH SHA256=${TESTDATA_SHA256}
  )
else()
# TESTDATA is set; use it.
  set(TESTDATA_DIR "${TESTDATA}")
  set(testdata_args )
endif()  

ADD_DEFINITIONS(-DTESTDATA_DIR=\"${TESTDATA_DIR}\")

# Create a target `testdata` to download and extract testdata as neccessary
externalproject_add(testdata
  PREFIX testdata
  EXCLUDE_FROM_ALL 1
  SOURCE_DIR "${TESTDATA_DIR}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""

  ${testdata_args}
)

###############################################################################
# Create a make_testdata_sum target
# This target will download testdata and print its checksum.  Useful for developers.
configure_file(
  "${CMAKE_SOURCE_DIR}/Modules/cmake_make_testdata_sum.cmake.in"
  "${CMAKE_BINARY_DIR}/cmake_make_testdata_sum.cmake"
  IMMEDIATE @ONLY)

add_custom_target(make_testdata_sum
  "${CMAKE_COMMAND}" -P "${CMAKE_BINARY_DIR}/cmake_make_testdata_sum.cmake"
  COMMENT "Calculating the checksum for testdata.tar.gz"
)

###############################################################################
# Determine the location of the executables/scripts and any data sets required for full package testing
set(DNG_DNM_EXE ${CMAKE_BINARY_DIR}/src/dng-dnm)
set(DNG_PHASER_EXE ${CMAKE_BINARY_DIR}/src/dng-phaser)
set(DNG_CALL_EXE ${CMAKE_BINARY_DIR}/src/dng-call)
set(DNG_LOGLIKE_EXE ${CMAKE_BINARY_DIR}/src/dng-loglike)
set(DNG_PILEUP_EXE ${CMAKE_BINARY_DIR}/src/dng-pileup)

###############################################################################
# Add tests for various build targets so we can depend on them
set(build_test_targets testdata dng-call dng-loglike dng-dnm dng-phaser dng-pileup)
foreach(target ${build_test_targets})
  add_test(Build.${target} "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${target})
endforeach()

###############################################################################
# AddCMakeTest
# A macro to create a test based on a .cmake.in file
# Example: AddCMakeTest(DngCall.Bam "" "Build.dng-call;Build.testdata")
#   - Creates a test called DngCall.Bam
#   - Makes that test depend on Build.dng-call and Build.testdata
#   - Processes the file ./DngCall/BamTest.cmake.in
#   - When the test is called, it calls the .cmake script with possible PREARGS and post ARGs

macro(AddCMakeTest NAME PREARGS DEPENDS)
  string(REPLACE "." "/" pathname "${NAME}")
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${pathname}Test.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${pathname}Test.cmake" @ONLY)
  add_test(NAME ${NAME}
    COMMAND ${CMAKE_COMMAND} ${PREARGS}
    -P "${CMAKE_CURRENT_BINARY_DIR}/${pathname}Test.cmake" ${ARGN})
  if("${DEPENDS}")
    set_tests_properties(${NAME} PROPERTIES DEPENDS "${DEPENDS}")
  endif()
endmacro()

###############################################################################
# Add Tests Here

add_test(Source.CheckBuildDir ${CMAKE_COMMAND} -E md5sum "${CMAKE_SOURCE_DIR}/build/.gitkeep")
set_tests_properties(Source.CheckBuildDir PROPERTIES
  PASS_REGULAR_EXPRESSION "d41d8cd98f00b204e9800998ecf8427e"
)

add_subdirectory(Unit)

AddCMakeTest(DngCall.Run "" "Build.dng-call")
AddCMakeTest(DngCall.Bam "" "Build.dng-call;Build.testdata")
AddCMakeTest(DngCall.Vcf "" "Build.dng-call;Build.testdata")
AddCmakeTest(DngCall.Cram "" "Build.dng-call;Build.testdata")
AddCmakeTest(DngCall.Tad "" "Build.dng-call;Build.testdata")

AddCMakeTest(DngLoglike.Run "" "Build.dng-loglike")

AddCMakeTest(DngDnm.Run "" "Build.dng-dnm")
AddCMakeTest(DngDnm.Auto "" "Build.dng-dnm;Build.testdata")

AddCMakeTest(DngPhaser.Run "" "Build.dng-phaser")
AddCMakeTest(DngPhaser.Data "" "Build.dng-phaser;Build.testdata")

AddCMakeTest(DngPileup.Run "" "Build.dng-pileup")
AddCMakeTest(DngPileup.Bam "" "Build.dng-pileup;Build.testdata")
