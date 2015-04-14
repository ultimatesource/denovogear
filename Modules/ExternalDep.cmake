INCLUDE(ExternalProject)

SET(EXT_PREFIX ext_deps)

SET(EXT_CFLAGS "${CMAKE_C_FLAGS}")
SET(EXT_LDFLAGS "${CMAKE_STATIC_LINKER_FLAGS}")
IF(CMAKE_BUILD_TYPE)
  STRING(TOUPPER "${CMAKE_BUILD_TYPE}" cmake_build_type_toupper)
  SET(EXT_CFLAGS "${EXT_CFLAGS} ${CMAKE_C_FLAGS_${cmake_build_type_toupper}}")
  SET(EXT_LDFLAGS "${EXT_LDFLAGS} ${CMAKE_STATIC_LINKER_FLAGS_${cmake_build_type_toupper}}")
ENDIF()

IF(NOT BUILD_EXTERNAL_PROJECTS)
  SET(REQ REQUIRED)
  SET(QUI "")
ELSE()
  SET(REQ QUIET)
  SET(QUI QUIET)
ENDIF()
SET(missing_ext_deps FALSE)

################################################################################
# BOOST
#

FIND_PACKAGE(Boost 1.47.0 REQUIRED COMPONENTS program_options filesystem system)

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
  	--with-graph
  	--disable-icu
  	--ignore-site-config
  	threading=multi
  	link=static
  	runtime-link=shared
  	variant=release
  	cxxflags=-std=c++11
  )

  ExternalProject_add(boost
    URL http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.bz2/download
    URL_MD5 1be49befbdd9a5ce9def2983ba3e7b76
    PREFIX ext_deps/boost
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ${boost_bootstrap}
    BUILD_COMMAND ${boost_build}
    INSTALL_COMMAND ""
  )
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

FIND_PACKAGE(Eigen3 ${QUI})

IF(BUILD_EXTERNAL_PROJECTS AND NOT EIGEN3_FOUND)
  ExternalProject_Add(eigen3
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
  SET(EIGEN3_EXT_TARGET eigen3)
  MESSAGE(STATUS "Building Eigen3 3.2.4 as external dependency")
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

FIND_PACKAGE(HTSLIB 1.2 ${REQ})

IF(BUILD_EXTERNAL_PROJECTS AND NOT HTSLIB_FOUND)
  ExternalProject_Add(htslib
    PREFIX "${EXT_PREFIX}/htslib"
    URL https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2
    URL_MD5 88eec909855abd98032bc2f9c3e83271
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND "./configure" "--prefix=<INSTALL_DIR>"
      "CC=${CMAKE_C_COMPILER}" "CFLAGS=${EXT_CFLAGS}" "LDFLAGS=${EXT_LDFLAGS}"
    BUILD_IN_SOURCE true
  )

  SET(HTSLIB_FOUND true)
  SET(HTSLIB_VERSION 1.2.1)
  SET(HTSLIB_LIBRARIES "${CMAKE_CURRENT_BINARY_DIR}/${EXT_PREFIX}/htslib/lib/libhts.a")
  SET(HTSLIB_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/${EXT_PREFIX}/htslib/include/")
  SET(HTSLIB_EXT_TARGET htslib)
  MESSAGE(STATUS "Building HTSLIB 1.2.1 as external dependency")
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
    "configuration variable\n  BUILD_EXTERNAL_PROJECTS=1\n"
    "E.g. cmake -DBUILD_EXTERNAL_PROJECTS=1 .\n")
ENDIF()