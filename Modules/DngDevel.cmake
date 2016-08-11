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
