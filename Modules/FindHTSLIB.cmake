#FIND_PACKAGE(PkgConfig)
#
#pkg_check_modules(HTSLIB htslib)

include(LibFindMacros)

libfind_pkg_detect(HTSLIB htslib FIND_PATH htslib/hts.h FIND_LIBRARY hts)

libfind_process(HTSLIB)

# set(_HTSLIB_SEARCHES)

# # Search HTSLIB_ROOT first if it is set.
# if(HTSLIB_ROOT)
#   set(_HTSLIB_SEARCH_ROOT PATHS ${HTSLIB_ROOT} NO_DEFAULT_PATH)
#   list(APPEND _HTSLIB_SEARCHES _HTSLIB_SEARCH_ROOT)
# endif()

# # Normal search.
# set(_HTSLIB_SEARCH_NORMAL)
# list(APPEND _HTSLIB_SEARCHES _HTSLIB_SEARCH_NORMAL)

# set(HTSLIB_NAMES hts)

# foreach(search ${_ZLIB_SEARCHES})
#   find_path(HTSLIB_INCLUDE_DIR NAMES htslib/hts.h        ${${search}} PATH_SUFFIXES include)
#   find_library(HTSLIB_LIBRARY  NAMES ${HTSLIB_NAMES} ${${search}} PATH_SUFFIXES lib)
# endforeach()

# mark_as_advanced(HTSLIB_LIBRARY HTSLIB_INCLUDE_DIR)

# # handle the QUIETLY and REQUIRED arguments and set HTSLIB_FOUND to TRUE if
# # all listed variables are TRUE
# include(FindPackageHandleStandardArgs)
# FIND_PACKAGE_HANDLE_STANDARD_ARGS(HTSLIB REQUIRED_VARS HTSLIB_LIBRARY HTSLIB_INCLUDE_DIR)

# if(HTSLIB_FOUND)
#     set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR})
#     set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY})
# endif()
