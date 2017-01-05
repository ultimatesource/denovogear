find_package(Rscript)

if (Rscript_FOUND)

  set(MUTMAP_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR}/mutmap)
  set(JS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/js)
  set(DIST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dist)
  set(BUNDLE_FILE ${CMAKE_CURRENT_BINARY_DIR}/bundle.js)
  
  # Note: order is important for dependencies here
  set(JS_FILES ${JS_DIR}/utils.js ${JS_DIR}/dng_vcf_parser.js
      ${JS_DIR}/ped_parser.js ${JS_DIR}/pedigr.js ${JS_DIR}/visuals.js
      ${JS_DIR}/dng_vis.js)
  
  add_custom_command(
      OUTPUT ${BUNDLE_FILE}
      COMMAND cat ${JS_FILES} > ${BUNDLE_FILE}
      DEPENDS ${JS_FILES}
      )
  
  add_custom_target(bundle ALL DEPENDS ${BUNDLE_FILE})
  
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/build.R"
    "${CMAKE_CURRENT_BINARY_DIR}/build.R" @ONLY)

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/layout_and_template.R"
    "${CMAKE_CURRENT_BINARY_DIR}/layout_and_template.R" @ONLY)

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/test_server.R"
    "${CMAKE_CURRENT_BINARY_DIR}/test_server.R" @ONLY)

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/autotest.sh"
    "${CMAKE_CURRENT_BINARY_DIR}/autotest.sh" @ONLY)
  
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/html/index.html"
    "${CMAKE_CURRENT_BINARY_DIR}/index.html" @ONLY)
  
  install (DIRECTORY DESTINATION ${MUTMAP_DIR})
  
  install (PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/build.R"
    DESTINATION ${CMAKE_INSTALL_LIBEXECDIR} RENAME dng-mutmap)

  install (FILES "${CMAKE_CURRENT_BINARY_DIR}/layout_and_template.R"
    DESTINATION ${MUTMAP_DIR})

  install (FILES ${BUNDLE_FILE} DESTINATION ${MUTMAP_DIR})
  
  install (FILES "${CMAKE_CURRENT_BINARY_DIR}/index.html"
    DESTINATION ${MUTMAP_DIR})

endif()
