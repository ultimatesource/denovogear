find_program(RSCRIPT_COMMAND Rscript DOC "Rscript executable.")

if (NOT RSCRIPT_COMMAND)
    message(STATUS "Rscript not found. It is required if you want to generate a visualization via dng mutmap. Please install R with Rscript")
endif()
