include(LibFindMacros)

IF(NOT ZLIB_FOUND)
  FIND_PACKAGE(ZLIB)
ENDIF()

libfind_pkg_detect(HTSLIB htslib FIND_PATH htslib/hts.h FIND_LIBRARY hts)

if(NOT HTSLIB_VERSION)
  if(HTSLIB_PKGCONF_VERSION)
    set(HTSLIB_VERSION "${HTSLIB_PKGCONF_VERSION}")
  else()
    FIND_PROGRAM(HTSLIB_EXECUTABLE NAMES htsfile)
    execute_process(
      COMMAND htsfile "--version"
      OUTPUT_VARIABLE hts_version_raw
    )
    string(REGEX MATCH "[0-9.]+" hts_version ${hts_version_raw})
    set(HTSLIB_VERSION ${hts_version})
  endif()
endif()

libfind_process(HTSLIB)


if(HTSLIB_FOUND)
  if(NOT TARGET HTSLIB::HTSLIB)
    add_library(HTSLIB::HTSLIB UNKNOWN IMPORTED)
    set_target_properties(HTSLIB::HTSLIB PROPERTIES
      IMPORTED_LOCATION "${HTSLIB_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES ZLIB::ZLIB
    )
    ADD_DEPENDENCIES(HTSLIB::HTSLIB ZLIB::ZLIB)
  endif()
endif(HTSLIB_FOUND)
