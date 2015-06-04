INCLUDE(ExternalProject)

SET(EXT_PREFIX ext_deps)

SET(EXT_CFLAGS "${CMAKE_C_FLAGS}")
SET(EXT_LDFLAGS "${CMAKE_STATIC_LINKER_FLAGS}")
IF(CMAKE_BUILD_TYPE)
  STRING(TOUPPER "${CMAKE_BUILD_TYPE}" cmake_build_type_toupper)
  SET(EXT_CFLAGS "${EXT_CFLAGS} ${CMAKE_C_FLAGS_${cmake_build_type_toupper}}")
  SET(EXT_LDFLAGS "${EXT_LDFLAGS} ${CMAKE_STATIC_LINKER_FLAGS_${cmake_build_type_toupper}}")

  ## Turn off debugging in libraries
  IF(cmake_build_type_toupper STREQUAL "RELEASE" OR
     cmake_build_type_toupper STREQUAL "RELWITHDEBINFO")
    ADD_DEFINITIONS("-DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG")
    SET(boost_variant variant=release)
  ELSE()
    SET(boost_variant variant=debug)
  ENDIF()
ENDIF()

IF(NOT BUILD_EXTERNAL_PROJECTS)
  SET(REQ REQUIRED)
  SET(QUI )
ELSE()
  SET(REQ )
  SET(QUI QUIET)
  STRING(TOUPPER "${BUILD_EXTERNAL_PROJECTS}" build_external_projects_toupper)
  IF(build_external_projects_toupper STREQUAL "FORCE")
    SET(BUILD_EXTERNAL_PROJECTS_FORCED ON)
  ENDIF()
ENDIF()
SET(missing_ext_deps FALSE)

IF(USE_STATIC_LIBS)
  SET(Boost_USE_STATIC_LIBS ON)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
ENDIF(USE_STATIC_LIBS)

################################################################################
# THREADS
#
SET(THREADS_PREFER_PTHREAD_FLAG ON)
FIND_PACKAGE(Threads)

################################################################################
# ZLIB
#

IF(NOT BUILD_EXTERNAL_PROJECTS_FORCED)
  FIND_PACKAGE(ZLIB)
ENDIF()

IF(BUILD_EXTERNAL_PROJECTS AND NOT ZLIB_FOUND)
  ExternalProject_add(ext_zlib
    URL http://zlib.net/zlib-1.2.8.tar.gz
    URL_MD5 44d667c142d7cda120332623eab69f40
    PREFIX ext_deps/zlib
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND "env" "CC=${CMAKE_C_COMPILER}" "CFLAGS=${EXT_CFLAGS}" "LDFLAGS=${EXT_LDFLAGS}"
      "./configure" "--prefix=<INSTALL_DIR>"
    BUILD_IN_SOURCE true
  )
  SET(ZLIB_FOUND TRUE)
  SET(ZLIB_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/${EXT_PREFIX}/zlib/include/")
  SET(ZLIB_LIBRARIES "${CMAKE_CURRENT_BINARY_DIR}/${EXT_PREFIX}/zlib/lib/libz.a")
  if(NOT TARGET ZLIB::ZLIB)
    add_library(ZLIB::ZLIB UNKNOWN IMPORTED)
    set_target_properties(ZLIB::ZLIB PROPERTIES
      IMPORTED_LOCATION "${ZLIB_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${ZLIB_INCLUDE_DIRS}")
    ADD_DEPENDENCIES(ZLIB::ZLIB ext_zlib)
    FILE(MAKE_DIRECTORY "${ZLIB_INCLUDE_DIRS}")
  endif()
  SET(ZLIB_VERSION_STRING "1.2.8")
  MESSAGE(STATUS "Building ZLIB ${ZLIB_VERSION_STRING} as external dependency")
ENDIF()

IF(ZLIB_FOUND)
ELSE()
  SET(missing_ext_deps TRUE)
ENDIF(ZLIB_FOUND)

################################################################################
# BOOST
#

IF(NOT BUILD_EXTERNAL_PROJECTS_FORCED)
  FIND_PACKAGE(Boost 1.47.0 ${REQ} COMPONENTS
    program_options
    filesystem
    system
    unit_test_framework
  )
ENDIF()

IF(BUILD_EXTERNAL_PROJECTS AND NOT Boost_FOUND)
  SET(boost_bootstrap "./bootstrap.sh")
  IF(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    LIST(APPEND boost_bootstrap --with-toolset=clang)
  ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    LIST(APPEND boost_bootstrap --with-toolset=gcc)
  endif()

  SET(boost_build "./b2" install
    --prefix=<INSTALL_DIR>
    --with-program_options
    --with-filesystem
    --with-system
    --with-test
    --disable-icu
    --ignore-site-config
    threading=multi
    link=static
    runtime-link=shared
    cxxflags=-std=c++11
    ${boost_variant}
  )

  ExternalProject_add(ext_boost
    URL http://downloads.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.tar.bz2
    URL_MD5 1be49befbdd9a5ce9def2983ba3e7b76
    PREFIX ext_deps/boost
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ${boost_bootstrap}
    BUILD_COMMAND ${boost_build}
    INSTALL_COMMAND ""
  )
  SET(boost_ext_libdir "${CMAKE_CURRENT_BINARY_DIR}/${EXT_PREFIX}/boost/lib")
  SET(Boost_FOUND TRUE)
  SET(Boost_VERSION 1.57)
  SET(Boost_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/${EXT_PREFIX}/boost/include/")

  SET(Boost_LIBRARIES)
  FOREACH(ext_boost_name PROGRAM_OPTIONS FILESYSTEM SYSTEM UNIT_TEST_FRAMEWORK)
    STRING(TOLOWER "${ext_boost_name}" ext_boost_lowname)
    SET(Boost_${ext_boost_name}_FOUND On)
    SET(Boost_${ext_boost_name}_LIBRARY "${boost_ext_libdir}/libboost_${ext_boost_lowname}.a")
    SET(Boost_LIBRARIES ${Boost_LIBRARIES} ${Boost_${ext_boost_name}_LIBRARY})
  ENDFOREACH()

  SET(Boost_LIBRARY_DIRS "")
  SET(BOOST_EXT_TARGET ext_boost)
  SET(Boost_USE_STATIC_LIBS TRUE)
  MESSAGE(STATUS "Building Boost ${Boost_VERSION} as external dependency")  
ENDIF()

IF(Boost_FOUND)
  ADD_DEFINITIONS(-DBOOST_ALL_NO_LIB -DBOOST_PROGRAM_OPTIONS_NO_LIB)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
  IF(NOT Boost_USE_STATIC_LIBS)
    ADD_DEFINITIONS(-DBOOST_DYN_LINK -DBOOST_PROGRAM_OPTIONS_DYN_LINK)
  ENDIF(NOT Boost_USE_STATIC_LIBS)
ELSE()
  SET(missing_ext_deps TRUE)
ENDIF(Boost_FOUND)

################################################################################
# EIGEN3
#

IF(NOT BUILD_EXTERNAL_PROJECTS_FORCED)
  FIND_PACKAGE(Eigen3 ${QUI})
ENDIF()

IF(BUILD_EXTERNAL_PROJECTS AND NOT EIGEN3_FOUND)
  ExternalProject_Add(ext_eigen3
    PREFIX "${EXT_PREFIX}/eigen3"
    URL http://bitbucket.org/eigen/eigen/get/3.2.4.tar.bz2
    URL_MD5 4c4b5ed9a388a1e475166d575af25477
    CMAKE_CACHE_ARGS "-DCMAKE_INSTALL_PREFIX:string=<INSTALL_DIR>"
     "-DCMAKE_CXX_COMPILER:filepath=${CMAKE_CXX_COMPILER}"
     "-DCMAKE_C_COMPILER:filepath=${CMAKE_C_COMPILER}"
  )
  SET(EIGEN3_FOUND TRUE)
  SET(EIGEN3_VERSION 3.2.4)
  SET(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/${EXT_PREFIX}/eigen3/include/eigen3/")
  SET(EIGEN3_EXT_TARGET ext_eigen3)
  MESSAGE(STATUS "Building Eigen3 ${EIGEN3_VERSION} as external dependency")
ENDIF()

IF(EIGEN3_FOUND)
  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR}/unsupported)
ELSE()
  SET(missing_ext_deps TRUE)
ENDIF(EIGEN3_FOUND)

################################################################################
# HTSLIB
#

IF(NOT BUILD_EXTERNAL_PROJECTS_FORCED)
  FIND_PACKAGE(HTSLIB 1.2 ${REQ} ${QUI})
ENDIF()

IF(BUILD_EXTERNAL_PROJECTS AND NOT HTSLIB_FOUND)
  GET_FILENAME_COMPONENT(zlib_lib_dir ${ZLIB_LIBRARIES} DIRECTORY)
  ExternalProject_Add(ext_htslib
    PREFIX "${EXT_PREFIX}/htslib"
    URL https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2
    URL_MD5 88eec909855abd98032bc2f9c3e83271
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND "./configure" "--prefix=<INSTALL_DIR>"
      "CC=${CMAKE_C_COMPILER}"
      "CFLAGS=${EXT_CFLAGS} -I${ZLIB_INCLUDE_DIRS}"
      "LDFLAGS=${EXT_LDFLAGS} -L${zlib_lib_dir}"
    BUILD_IN_SOURCE true
  )

  SET(HTSLIB_FOUND true)
  SET(HTSLIB_VERSION 1.2.1)
  SET(HTSLIB_LIBRARIES "${CMAKE_CURRENT_BINARY_DIR}/${EXT_PREFIX}/htslib/lib/libhts.a")
  SET(HTSLIB_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/${EXT_PREFIX}/htslib/include/")
  SET(HTSLIB_EXT_TARGET ext_htslib)
  IF(ZLIB_FOUND)
    ADD_DEPENDENCIES(ext_htslib, ZLIB::ZLIB)
  ENDIF(ZLIB_FOUND)
  MESSAGE(STATUS "Building HTSLIB ${HTSLIB_VERSION} as external dependency")
ENDIF()

IF(HTSLIB_FOUND)
  INCLUDE_DIRECTORIES(${HTSLIB_INCLUDE_DIRS})
ELSE()
  SET(missing_ext_deps TRUE)
ENDIF(HTSLIB_FOUND)

################################################################################
# Help Message
#

IF(NOT BUILD_EXTERNAL_PROJECTS AND missing_ext_deps)
  MESSAGE(SEND_ERROR "ERROR: Some dependecies could not be found on your system. "
    "To download and build the missing dependecies please set the following "
    "configuration variable\n    BUILD_EXTERNAL_PROJECTS=1\n"
    "E.g. cmake -DBUILD_EXTERNAL_PROJECTS=1 .\n")
ENDIF()
