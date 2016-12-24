find_package(Rscript)

set(MUTMAP_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR}/mutmap)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/mutmap.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mutmap.sh" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/dist/build.R"
  "${CMAKE_CURRENT_BINARY_DIR}/build.R" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/dist/bundle.js"
  "${CMAKE_CURRENT_BINARY_DIR}/bundle.js" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/dist/index.html"
  "${CMAKE_CURRENT_BINARY_DIR}/index.html" @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/dist/pedigree_and_layout.R"
  "${CMAKE_CURRENT_BINARY_DIR}/pedigree_and_layout.R" @ONLY)

install (DIRECTORY DESTINATION ${MUTMAP_DIR})

install (PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/mutmap.sh"
  DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} RENAME dng-mutmap)

install (FILES "${CMAKE_CURRENT_BINARY_DIR}/build.R"
  DESTINATION ${MUTMAP_DIR})

install (FILES "${CMAKE_CURRENT_BINARY_DIR}/bundle.js"
  DESTINATION ${MUTMAP_DIR})

install (FILES "${CMAKE_CURRENT_BINARY_DIR}/index.html"
  DESTINATION ${MUTMAP_DIR})

install (FILES "${CMAKE_CURRENT_BINARY_DIR}/pedigree_and_layout.R"
  DESTINATION ${MUTMAP_DIR})

