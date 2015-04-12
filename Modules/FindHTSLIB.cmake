include(LibFindMacros)

libfind_pkg_detect(HTSLIB htslib FIND_PATH htslib/hts.h FIND_LIBRARY hts)

if(NOT HTSLIB_VERSION)
	if(HTSLIB_PKGCONF_VERSION)
		set(HTSLIB_VERSION "${HTSLIB_PKGCONF_VERSION}")
	endif()
endif()

libfind_process(HTSLIB)
