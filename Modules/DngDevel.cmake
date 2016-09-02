# This CMake File defines several useful developer options

SET(DNG_DEVEL_ENABLE_GPERFTOOLS OFF CACHE BOOL "Enable profiling with gperftools.")

SET(dng_devel_LIBRARIES)

if(DNG_DEVEL_ENABLE_GPERFTOOLS)
  find_package(Gperftools COMPONENTS profiler)
  if(GPERFTOOLS_FOUND)
    message(STATUS "DNG_DEVEL: Profiling with gperftools enabled. Use CPUPROFILE environmental variable to turn on profiling and specify output file.")
    set(dng_devel_LIBRARIES ${dng_devel_LIBRARIES} GPERFTOOLS::GPERFTOOLS)
  else()
    message(FATAL_ERROR "Gperftools was not found. Please disable the flag DNG_DEVEL_ENABLE_GPERFTOOLS and try again.")
  endif()
endif()

SET(DNG_DEVEL_ENABLE_COVERAGE_REPORT OFF CACHE BOOL "Enable code coverage reporting.")

if (DNG_DEVEL_ENABLE_COVERAGE_REPORT)
  ## Only compatible with debug builds
  if(CMAKE_BUILD_TYPE)
    string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
    if(NOT cmake_build_type_tolower STREQUAL "debug")
        message(FATAL_ERROR "Unsupported build type \"${CMAKE_BUILD_TYPE}\". DNG_DEVEL_ENABLE_COVERAGE_REPORT can only be used with a debug build")
    else()
        SET(COVERAGE_FLAGS --coverage)
        SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}")
        SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_FLAGS}")
    endif()
  endif()
endif()
