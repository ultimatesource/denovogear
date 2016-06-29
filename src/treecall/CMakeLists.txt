find_package(PythonInterp 2.7)
find_package(PythonLibs 2.7)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/treecall.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/treecall.py" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/treecall/setup.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/treecall/setup.py" @ONLY)

set(TCOUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build/treecall-timestamp")

add_custom_command(OUTPUT ${TCOUTPUT}
  COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/treecall/setup.py" build_py -q
  COMMAND ${CMAKE_COMMAND} -E touch ${TCOUTPUT}
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/treecall/__init__.py"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(treecall ALL DEPENDS ${TCOUTPUT})

install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/treecall.py"
  DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} RENAME dng-treecall)

install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/treecall/setup.py install_lib -q --install-dir \${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
