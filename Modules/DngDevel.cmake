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
        message(FATAL_ERROR "Unsupported build type \"${CMAKE_BUILD_TYPE}\". DNG_DEVEL_ENABLE_COVERAGE_REPORT can only be used with a debug build.")
    else()
        message(STATUS "DNG_DEVEL: Coverage report enabled.")
        SET(COVERAGE_FLAGS --coverage)
        SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${COVERAGE_FLAGS}")
        SET(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${COVERAGE_FLAGS}")
        SET(CMAKE_EXE_LINKER_FLAGS_DEBUG  "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${COVERAGE_FLAGS}")
        SET(CMAKE_MODULE_LINKER_FLAGS_DEBUG  "${CMAKE_MODULE_LINKER_FLAGS_DEBUG} ${COVERAGE_FLAGS}")
        SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG  "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} ${COVERAGE_FLAGS}")
    endif()
  endif()
endif()
