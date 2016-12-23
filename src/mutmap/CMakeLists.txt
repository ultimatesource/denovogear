find_package(PythonLibs 2.7)
find_package(Rscript)

set(MUTMAP_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR}/mutmap)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/mutmap.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mutmap.py" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pedigree_and_layout.R"
  "${CMAKE_CURRENT_BINARY_DIR}/pedigree_and_layout.R" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/index.html"
  "${CMAKE_CURRENT_BINARY_DIR}/index.html" @ONLY)

install (DIRECTORY DESTINATION ${MUTMAP_DIR})

install (PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/mutmap.py"
  DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} RENAME dng-mutmap)

install (FILES "${CMAKE_CURRENT_BINARY_DIR}/pedigree_and_layout.R"
  DESTINATION ${MUTMAP_DIR})

install (FILES "${CMAKE_CURRENT_BINARY_DIR}/index.html"
  DESTINATION ${MUTMAP_DIR})
