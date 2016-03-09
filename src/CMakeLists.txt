INCLUDE(CheckIncludeFile)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)

INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")
INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/include")
INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/dnm")

ADD_EXECUTABLE(dng-call
  dng-call.cc
  lib/call.cc lib/likelihood.cc lib/newick.cc lib/pedigree.cc lib/peeling.cc
  lib/mutation.cc lib/stats.cc lib/regions.cc
  lib/mutation_stats.cc lib/find_mutation.cc
)

# Some compilers have incomplete c++11 support provide a list of
# common features to toss an error if someone tries to use
# and early compiler
TARGET_COMPILE_FEATURES(dng-call PRIVATE
  cxx_auto_type
  cxx_constexpr
  cxx_decltype
  cxx_deleted_functions
  cxx_lambdas
  cxx_nonstatic_member_init
  cxx_nullptr
  cxx_range_for
  cxx_rvalue_references
  cxx_strong_enums
  cxx_trailing_return_types
  cxx_attributes
  cxx_reference_qualified_functions
  cxx_decltype_incomplete_return_types
)

ADD_EXECUTABLE(dng-dnm
  dng-dnm.cc
  dnm/snpLike.cc dnm/indelLike.cc dnm/pairLike.cc dnm/makeLookup.cc
  dnm/pedParser.cc dnm/bcf2Qcall.cc dnm/bcf2Paired.cc
)

ADD_EXECUTABLE(dng-phaser
  dng-phaser.cc
)

CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/dng.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/dng.sh"
  @ONLY)

CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/version.h"
  @ONLY)

find_package(PythonInterp 2.7 REQUIRED)
find_package(PythonLibs 2.7 REQUIRED)
  
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/treecall2.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/treecall2.py" @ONLY)
  
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/treecall/setup.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/treecall/setup.py" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/py-editdist/setup.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/py-editdist/setup.py" @ONLY)

add_custom_target(treecall ALL ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/treecall/setup.py build) 

add_custom_target(py-ed ALL ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/py-editdist/setup.py build)

IF(USE_STATIC_LIBSTDCXX)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++")
ENDIF()

IF(USE_STATIC_FLAG)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  SET_TARGET_PROPERTIES(dng-dnm    PROPERTIES LINK_SEARCH_START_STATIC On)
  SET_TARGET_PROPERTIES(dng-phaser PROPERTIES LINK_SEARCH_START_STATIC On)
  SET_TARGET_PROPERTIES(dng-call   PROPERTIES LINK_SEARCH_START_STATIC On)

  SET_TARGET_PROPERTIES(dng-dnm    PROPERTIES LINK_SEARCH_END_STATIC On)
  SET_TARGET_PROPERTIES(dng-phaser PROPERTIES LINK_SEARCH_END_STATIC On)
  SET_TARGET_PROPERTIES(dng-call   PROPERTIES LINK_SEARCH_END_STATIC On)
ENDIF(USE_STATIC_FLAG)

TARGET_LINK_LIBRARIES(dng-dnm Boost::PROGRAM_OPTIONS Boost::FILESYSTEM Boost::SYSTEM)
TARGET_LINK_LIBRARIES(dng-dnm HTSLIB::HTSLIB Threads::Threads)
TARGET_LINK_LIBRARIES(dng-dnm EIGEN3::EIGEN3)

TARGET_LINK_LIBRARIES(dng-phaser Boost::PROGRAM_OPTIONS Boost::FILESYSTEM Boost::SYSTEM)
TARGET_LINK_LIBRARIES(dng-phaser HTSLIB::HTSLIB Threads::Threads)

TARGET_LINK_LIBRARIES(dng-call HTSLIB::HTSLIB Threads::Threads)
TARGET_LINK_LIBRARIES(dng-call EIGEN3::EIGEN3)
TARGET_LINK_LIBRARIES(dng-call Boost::PROGRAM_OPTIONS Boost::FILESYSTEM Boost::SYSTEM)

IF(DEVEL_MODE)
  TARGET_LINK_LIBRARIES(dng-call Boost::TIMER)
ENDIF()

IF(TARGET ext_eigen3)
  ADD_DEPENDENCIES(dng-call ext_eigen3)
  ADD_DEPENDENCIES(dng-dnm ext_eigen3)
ENDIF()

################################################################################
# Install

INSTALL(TARGETS dng-dnm RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} ${PROJ_INSTALL_EXTRA_ARGS})
INSTALL(TARGETS dng-phaser RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} ${PROJ_INSTALL_EXTRA_ARGS})
INSTALL(TARGETS dng-call RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} ${PROJ_INSTALL_EXTRA_ARGS})
INSTALL(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/dng.sh"
  DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME dng)
INSTALL(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/treecall2.py"
  DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} RENAME dng-treecall)
install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/treecall/setup.py install --install-purelib \${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR})")
install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/py-editdist/setup.py install --install-lib \${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR})")
