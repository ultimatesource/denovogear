include(CheckIncludeFile)
include(CheckFunctionExists)
include(CheckLibraryExists)

add_subdirectory(lib)
add_subdirectory(treecall)

include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/dng.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/dng.sh"
  @ONLY)
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/dng.sh" DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME dng)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/version.h"
  @ONLY)

if(USE_STATIC_LIBSTDCXX)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++")
endif()

if(USE_STATIC_FLAG)
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

# This is a list of binaries that will be compiled to run with dng
set(MODERN_DNG call loglike phaser dnm pileup) 

foreach(task IN LISTS MODERN_DNG)
  add_executable(dng-${task} dng-${task}.cc)
  if(USE_STATIC_FLAG)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    set_target_properties(dng-${task} PROPERTIES LINK_SEARCH_START_STATIC On)
    set_target_properties(dng-${task} PROPERTIES LINK_SEARCH_END_STATIC On)
  endif()
  target_link_libraries(dng-${task} libdng)
  target_link_libraries(dng-${task} ${dng_devel_LIBRARIES})

  install(TARGETS dng-${task} RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} ${PROJ_INSTALL_EXTRA_ARGS})
endforeach()

if(DEVEL_MODE)
  target_link_libraries(dng-call Boost::TIMER)
endif()

# Add extra configuration options for DNM
target_include_directories(dng-dnm BEFORE PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/dnm")

target_sources(dng-dnm PRIVATE
  dnm/snpLike.cc dnm/indelLike.cc dnm/pairLike.cc dnm/makeLookup.cc
  dnm/pedParser.cc dnm/bcf2Qcall.cc dnm/bcf2Paired.cc)

