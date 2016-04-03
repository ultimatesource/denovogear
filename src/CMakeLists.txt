include(CheckIncludeFile)
include(CheckFunctionExists)
include(CheckLibraryExists)

add_subdirectory(lib)
add_subdirectory(treecall)

include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(dng-call dng-call.cc)

add_executable(dng-dnm dng-dnm.cc
  dnm/snpLike.cc dnm/indelLike.cc dnm/pairLike.cc dnm/makeLookup.cc
  dnm/pedParser.cc dnm/bcf2Qcall.cc dnm/bcf2Paired.cc
)

target_include_directories(dng-dnm BEFORE PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/dnm")

add_executable(dng-phaser dng-phaser.cc)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/dng.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/dng.sh"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/version.h"
  @ONLY)
  
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/py-editdist/setup.py.in"
#  "${CMAKE_CURRENT_BINARY_DIR}/py-editdist/setup.py" @ONLY)

#add_custom_target(py-editdist ALL ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/py-editdist/setup.py build)

if(USE_STATIC_LIBSTDCXX)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++")
endif()

if(USE_STATIC_FLAG)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  set_target_properties(dng-dnm    PROPERTIES LINK_SEARCH_START_STATIC On)
  set_target_properties(dng-phaser PROPERTIES LINK_SEARCH_START_STATIC On)
  set_target_properties(dng-call   PROPERTIES LINK_SEARCH_START_STATIC On)

  set_target_properties(dng-dnm    PROPERTIES LINK_SEARCH_END_STATIC On)
  set_target_properties(dng-phaser PROPERTIES LINK_SEARCH_END_STATIC On)
  set_target_properties(dng-call   PROPERTIES LINK_SEARCH_END_STATIC On)
endif(USE_STATIC_FLAG)

target_link_libraries(dng-call libdng)
target_link_libraries(dng-dnm libdng)
target_link_libraries(dng-phaser libdng)

if(DEVEL_MODE)
  target_link_libraries(dng-call Boost::TIMER)
endif()

################################################################################
# Install

install(TARGETS dng-dnm RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} ${PROJ_INSTALL_EXTRA_ARGS})
install(TARGETS dng-phaser RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} ${PROJ_INSTALL_EXTRA_ARGS})
install(TARGETS dng-call RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} ${PROJ_INSTALL_EXTRA_ARGS})
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/dng.sh" DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME dng)

#install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/py-editdist/setup.py install --install-lib \${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR})")
