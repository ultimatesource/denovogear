################################################################################
#  Denovogear Project Configuration
#
## Project DeNovoGear
CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0 FATAL_ERROR)
CMAKE_POLICY(SET CMP0025 NEW)
# CMAKE_POLICY(SET CMP0022 NEW)
# CMAKE_POLICY(SET CMP0079 NEW)

## Important Variables
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Modules")
SET(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")
SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

## Compiling options
SET(USE_STATIC_LIBS OFF CACHE BOOL "Link against static libraries.")
SET(USE_STATIC_FLAG OFF CACHE BOOL "Use -static flag to create completely static binaries.")
SET(BUILD_EXTERNAL_PROJECTS OFF CACHE BOOL "Build external dependencies if they are not found.")
SET(DEVEL_MODE OFF CACHE BOOL "Enable developer mode.")
SET(EXTRA_DEFS CACHE STRING "Extra Compiler Definations")

################################################################################
# Init Project

PROJECT(DeNovoGear CXX C)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
#SET(CMAKE_CXX_EXTENSIONS ON)

# AppleClang doesn't know how to do these things until CMake 3.2
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" AND
  "${CMAKE_VERSION}" VERSION_LESS 3.2)
  SET(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
  SET(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++11")
ENDIF()

SET(cxx_std_flag "")
IF(CMAKE_CXX_STANDARD_REQUIRED)
  IF(CMAKE_CXX_EXTENSIONS)
    SET(cxx_std_flag ${CMAKE_CXX${CMAKE_CXX_STANDARD}_EXTENSION_COMPILE_OPTION} )
  ELSE()
    SET(cxx_std_flag ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} )
  ENDIF()
ENDIF()

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

## Guard against bad build-type strings
IF(CMAKE_BUILD_TYPE)
  string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
  if(    NOT cmake_build_type_tolower STREQUAL "debug"
     AND NOT cmake_build_type_tolower STREQUAL "release"
     AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo"
     AND NOT cmake_build_type_tolower STREQUAL "none")
    message(FATAL_ERROR "Unsupported build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo, and None (case-insensitive).")
  endif()
  message(STATUS "Compiling with build type '${CMAKE_BUILD_TYPE}'.")
  if(cmake_build_type_tolower STREQUAL "debug")
    if(NOT DEFINED CMAKE_EXPORT_COMPILE_COMMANDS)
      set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    endif()
  endif()
ENDIF()


################################################################################
# Always make a full statement RPATH
# See http://www.itk.org/Wiki/CMake_RPATH_handling

SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}" is_system_dir)
IF("${is_system_dir}" STREQUAL "-1")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
ENDIF("${is_system_dir}" STREQUAL "-1")
################################################################################

IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
  ADD_DEFINITIONS(-D_GLIBCXX_USE_C99 -D_GLIBCXX_USE_C99_MATH -D_GLIBCXX_USE_C99_MATH_TR1 -D_DECLARE_C99_LDBL_MATH)
ENDIF()

## Specify Project Information
SET(DNG_NAME "DeNovoGear")
SET(DNG_NAME_FOR_PKG_FILE "denovogear")
SET(DNG_URL "https://github.com/denovogear/denovogear")

IF(DEVEL_MODE)
  ADD_DEFINITIONS(-DDNG_DEVEL)
ENDIF()

IF(EXTRA_DEFS)
  SEPARATE_ARGUMENTS(deflist UNIX_COMMAND "${EXTRA_DEFS}")
  FOREACH(def IN LISTS deflist)
    ADD_DEFINITIONS(-D${def})
  ENDFOREACH()
ENDIF()

## Turn on unit and executable testing
ENABLE_TESTING()

## Enable Developer CMake Flags
INCLUDE(DngDevel)

################################################################################
# Version Information

## Detect Project Version Information
FIND_PACKAGE(Git)
SET(GIT_COMMIT_ID "unknown")
IF(GIT_FOUND)
  EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res_var OUTPUT_VARIABLE GIT_COM_ID
     ERROR_QUIET)
  IF(${res_var} EQUAL 0 )
   	STRING( REPLACE "\n" "" GIT_COMMIT_ID "${GIT_COM_ID}" )
  ENDIF()
ENDIF()

# try directory name
IF(NOT GIT_COMMIT_ID AND PROJECT_SOURCE_DIR MATCHES "${DNG_NAME_FOR_PKG_FILE}-([-a-zA-Z0-9.]*)$" )
  SET(GIT_COMMIT_ID "${CMAKE_MATCH_1}")
ENDIF()

# extract project info
IF(GIT_COMMIT_ID MATCHES
  "^v?([0-9]+)\\.([0-9]+)\\.?([0-9]+)?((-[a-zA-Z]+)?(-[0-9]+)?(-g[0-9A-Fa-f]+)?(-dirty)?)$" )
  IF(NOT DEFINED DNG_VERSION_MAJOR)
    SET(DNG_VERSION_MAJOR "${CMAKE_MATCH_1}")
  ENDIF()
  IF(NOT DEFINED DNG_VERSION_MINOR)
    SET(DNG_VERSION_MINOR "${CMAKE_MATCH_2}")
  ENDIF()
  IF(NOT DEFINED DNG_VERSION_PATCH)
    SET(DNG_VERSION_PATCH "${CMAKE_MATCH_3}")
  ENDIF()
  IF(NOT DEFINED DNG_VERSION_COMMIT)
    STRING(REGEX REPLACE "^-" "" DNG_VERSION_COMMIT "${CMAKE_MATCH_4}")
  ENDIF()
ELSE()
  # We were unable to determine project information, use the default values
  INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in")
ENDIF()

SET(DNG_VERSION "${DNG_VERSION_MAJOR}.${DNG_VERSION_MINOR}")
SET(DNG_VERSION_SHORT "${DNG_VERSION_MAJOR}.${DNG_VERSION_MINOR}")
IF(DEFINED DNG_VERSION_PATCH)
  SET(DNG_VERSION "${DNG_VERSION}.${DNG_VERSION_PATCH}")
  SET(DNG_VERSION_SHORT "${DNG_VERSION_SHORT}.${DNG_VERSION_PATCH}")
ENDIF()
IF(DNG_VERSION_COMMIT)
  SET(DNG_VERSION "${DNG_VERSION}-${DNG_VERSION_COMMIT}")
  SET(DNG_VERSION_SHORT "${DNG_VERSION_SHORT}-${DNG_VERSION_COMMIT}")	
ENDIF()

MESSAGE(STATUS "${DNG_NAME} version: ${DNG_VERSION}.")

configure_file(DeNovoGear.sublime-project DeNovoGear.sublime-project COPYONLY)

################################################################################
# Packaging Information

SET(CMAKE_DATA_DIR "share/${DNG_NAME_FOR_PKG_FILE}" CACHE STRING
  "Install location for data (relative to prefix).")
SET(CMAKE_DOC_DIR "doc/${DNG_NAME_FOR_PKG_FILE}" CACHE STRING
  "Install location for documentation (relative to prefix).")
SET(CMAKE_MAN_DIR "man" CACHE STRING
  "Install location for man pages (relative to prefix).")
MARK_AS_ADVANCED(CMAKE_DATA_DIR CMAKE_DOC_DIR CMAKE_MAN_DIR)

SET(CPACK_SOURCE_IGNORE_FILES
  "/CVS/" "/\\\\.svn/"  "/\\\\.bzr/" "/\\\\.hg/" "/\\\\.git/" "\\\\.swp$"
  "\\\\.#"  "/#"  ".*~$"
  "/build/[^.]"
  "/CMakeFiles/"  "CMakeCache\\\\.txt$"
  "CPack.*Config\\\\.cmake$"  "cmake_install\\\\.cmake$"
  "install_manifest\\\\.txt$"
  "_CPACK_PACKAGES"  "_CPack_Packages"
  "\\\\.vcproj$"  "\\\\.dir$"  "\\\\.ncb$"  "\\\\.sln$"  "\\\\.suo$"
  "Makefile$" "\\\\.ilk$"  "\\\\.pdb$"
  "\\\\.a$"
)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Discover de-novo mutations from deep-sequencing of related individuals.")
SET(CPACK_PACKAGE_VENDOR "Conrad and Cartwright Labs")
SET(CPACK_PACKAGE_EXECUTABLES "denovogear" "DeNovoGear")

SET(CPACK_PACKAGE_INSTALL_DIRECTORY ${DNG_NAME_PKG})
SET(CPACK_PACKAGE_NAME ${DNG_NAME})
SET(CPACK_PACKAGE_NAME_FILE ${DNG_NAME_FOR_PKG_FILE})
SET(CPACK_PACKAGE_VERSION ${DNG_VERSION})
SET(CPACK_PACKAGE_VERSION_MAJOR ${DNG_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${DNG_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${DNG_VERSION_PATCH})
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING.txt")

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME_FILE}-v${CPACK_PACKAGE_VERSION}")

IF(NOT DEFINED CPACK_SYSTEM_NAME)
  SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
ENDIF(NOT DEFINED CPACK_SYSTEM_NAME)

IF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
  SET(CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
ENDIF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)

SET(CPACK_GENERATOR "TGZ;TBZ2")
SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2") 
SET(CPACK_STRIP_FILES "libexec/dng/*")

INCLUDE(GNUInstallDirs)
SET(CMAKE_INSTALL_LIBEXECDIR  "${CMAKE_INSTALL_LIBEXECDIR}/${DNG_NAME_FOR_PKG_FILE}")

INCLUDE(CPack)

INSTALL(FILES README.md COPYING.txt DESTINATION ${CMAKE_DOC_DIR})

INCLUDE(ExternalProject)
#INCLUDE(ExternalDep)

SET(THREADS_PREFER_PTHREAD_FLAG ON)
FIND_PACKAGE(Threads)
FIND_PACKAGE(ZLIB)
FIND_PACKAGE(Boost 1.47.0 REQUIRED COMPONENTS
  program_options
  filesystem
  system
  unit_test_framework
)
FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(HTSLIB REQUIRED)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(Tests)
ADD_SUBDIRECTORY(doc)


########### Add uninstall target ###############
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/Modules/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake") 

#### NOTES ####

## Build for 10.7+ on OSX
# CXXFLAGS=-stdlib=libc++ cmake -DCMAKE_BUILD_TYPE=Release \
#   -DUSE_STATIC_LIBS=on -DCMAKE_OSX_DEPLOYMENT_TARGET=10.7 ..
