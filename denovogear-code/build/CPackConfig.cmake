# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. Example variables are:
#   CPACK_GENERATOR                     - Generator used to create package
#   CPACK_INSTALL_CMAKE_PROJECTS        - For each project (path, name, component)
#   CPACK_CMAKE_GENERATOR               - CMake Generator used for the projects
#   CPACK_INSTALL_COMMANDS              - Extra commands to install components
#   CPACK_INSTALLED_DIRECTORIES           - Extra directories to install
#   CPACK_PACKAGE_DESCRIPTION_FILE      - Description file for the package
#   CPACK_PACKAGE_DESCRIPTION_SUMMARY   - Summary of the package
#   CPACK_PACKAGE_EXECUTABLES           - List of pairs of executables and labels
#   CPACK_PACKAGE_FILE_NAME             - Name of the package generated
#   CPACK_PACKAGE_ICON                  - Icon used for the package
#   CPACK_PACKAGE_INSTALL_DIRECTORY     - Name of directory for the installer
#   CPACK_PACKAGE_NAME                  - Package project name
#   CPACK_PACKAGE_VENDOR                - Package project vendor
#   CPACK_PACKAGE_VERSION               - Package project version
#   CPACK_PACKAGE_VERSION_MAJOR         - Package project version (major)
#   CPACK_PACKAGE_VERSION_MINOR         - Package project version (minor)
#   CPACK_PACKAGE_VERSION_PATCH         - Package project version (patch)

# There are certain generator specific ones

# NSIS Generator:
#   CPACK_PACKAGE_INSTALL_REGISTRY_KEY  - Name of the registry key for the installer
#   CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS - Extra commands used during uninstall
#   CPACK_NSIS_EXTRA_INSTALL_COMMANDS   - Extra commands used during install


SET(CPACK_BINARY_BUNDLE "")
SET(CPACK_BINARY_CYGWIN "")
SET(CPACK_BINARY_DEB "")
SET(CPACK_BINARY_DRAGNDROP "")
SET(CPACK_BINARY_NSIS "")
SET(CPACK_BINARY_OSXX11 "")
SET(CPACK_BINARY_PACKAGEMAKER "")
SET(CPACK_BINARY_RPM "")
SET(CPACK_BINARY_STGZ "")
SET(CPACK_BINARY_TBZ2 "")
SET(CPACK_BINARY_TGZ "")
SET(CPACK_BINARY_TZ "")
SET(CPACK_BINARY_ZIP "")
SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")
SET(CPACK_COMPONENTS_ALL "")
SET(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
SET(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
SET(CPACK_CREATE_DESKTOP_LINKS "denovogear")
SET(CPACK_GENERATOR "PackageMaker")
SET(CPACK_INSTALL_CMAKE_PROJECTS "/Users/dclab/conradlab/denovogear/denovogear-code/denovogear-code/build;DeNovoGear;ALL;/")
SET(CPACK_INSTALL_PREFIX "/usr/local/")
SET(CPACK_MODULE_PATH "/Users/dclab/conradlab/denovogear/denovogear-code/denovogear-code/Modules")
SET(CPACK_NSIS_DISPLAY_NAME "DeNovoGear")
SET(CPACK_NSIS_INSTALLER_ICON_CODE "")
SET(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
SET(CPACK_NSIS_PACKAGE_NAME "DeNovoGear")
SET(CPACK_OUTPUT_CONFIG_FILE "/Users/dclab/conradlab/denovogear/denovogear-code/denovogear-code/build/CPackConfig.cmake")
SET(CPACK_PACKAGE_DEFAULT_LOCATION "/usr/local/")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "/Users/dclab/conradlab/denovogear/denovogear-code/denovogear-code/README.txt")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Discover de novo mutations from deep-sequencing of related individuals.")
SET(CPACK_PACKAGE_EXECUTABLES "denovogear;DeNovoGear")
SET(CPACK_PACKAGE_FILE_NAME "denovogear-0.5-1-g8243639-Darwin-i386")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "DeNovoGear")
SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "DeNovoGear")
SET(CPACK_PACKAGE_NAME "DeNovoGear")
SET(CPACK_PACKAGE_NAME_FILE "denovogear")
SET(CPACK_PACKAGE_RELOCATABLE "true")
SET(CPACK_PACKAGE_VENDOR "Conrad and Cartwright Labs")
SET(CPACK_PACKAGE_VERSION "0.5-1-g8243639")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "5")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PACKAGING_INSTALL_PREFIX "/denovogear 0.5-1-g8243639")
SET(CPACK_RESOURCE_FILE_LICENSE "/Users/dclab/conradlab/denovogear/denovogear-code/denovogear-code/COPYING.txt")
SET(CPACK_RESOURCE_FILE_README "/Applications/CMake 2.8-7.app/Contents/share/cmake-2.8/Templates/CPack.GenericDescription.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "/Applications/CMake 2.8-7.app/Contents/share/cmake-2.8/Templates/CPack.GenericWelcome.txt")
SET(CPACK_SET_DESTDIR "OFF")
SET(CPACK_SOURCE_CYGWIN "")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~$;/\\.git/;/CMakeFiles/;CMakeCache\\.txt$;CPack.*Config\\.cmake$;cmake_install\\.cmake$;install_manifest\\.txt$;_CPACK_PACKAGES;_CPack_Packages;\\.vcproj$;\\.dir$;\\.ncb$;\\.sln$;\\.suo$;Makefile$;\\.ilk$;\\.pdb$;/releng/;\\.a$;denovogear[-]0;denovogear$;denovogear\\.exe")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/Users/dclab/conradlab/denovogear/denovogear-code/denovogear-code/build/CPackSourceConfig.cmake")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "denovogear-0.5-1-g8243639")
SET(CPACK_SOURCE_TBZ2 "")
SET(CPACK_SOURCE_TGZ "")
SET(CPACK_SOURCE_TZ "")
SET(CPACK_SOURCE_ZIP "")
SET(CPACK_SYSTEM_NAME "Darwin-i386")
SET(CPACK_TOPLEVEL_TAG "Darwin-i386")
