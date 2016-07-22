find_package(PythonInterp 2.7)

include(ExternalProject)

set(PREFIX_DIR "${CMAKE_CURRENT_BINARY_DIR}/ext_deps/python")
set(PKG_SUBDIR "lib/python2.7/site-packages")
set(BIN_SUBDIR "bin")
set(BIN_DIR "${PREFIX_DIR}/${BIN_SUBDIR}")
set(PKG_DIR "${PREFIX_DIR}/${PKG_SUBDIR}")
set(PKG_DEST ${CMAKE_INSTALL_PREFIX}/python/${PKG_SUBDIR})
set(BIN_DEST ${CMAKE_INSTALL_PREFIX}/python/${BIN_SUBDIR})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/treecall.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/treecall.sh" @ONLY)

file(MAKE_DIRECTORY ${PKG_DIR})

ExternalProject_add(
  ext_treecall
  URL ${CMAKE_CURRENT_SOURCE_DIR}/treecall
  PREFIX ${PREFIX_DIR}
  PATCH_COMMAND ""
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE true
  BUILD_COMMAND ""
  INSTALL_COMMAND PYTHONPATH=${PKG_DIR} ${PYTHON_EXECUTABLE} setup.py install --prefix ${PREFIX_DIR}
  DEPENDS ${dependencies}
  )


install(
  DIRECTORY ${PKG_DIR}/
  DESTINATION ${PKG_DEST}/
  USE_SOURCE_PERMISSIONS
  )

install(
  DIRECTORY ${BIN_DIR}/
  DESTINATION ${BIN_DEST}/
  USE_SOURCE_PERMISSIONS
  )

install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/treecall.sh"
  DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} RENAME dng-treecall)
