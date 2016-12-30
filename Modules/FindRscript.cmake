find_program(RSCRIPT_COMMAND Rscript DOC "Rscript executable.")

#find_package_handle_standard_args(Rscript DEFAULT_MSG RSCRIPT_COMMAND)
find_package_handle_standard_args(
  Rscript
  "Could not find Rscript, which is required for mutmap visualization" 
  RSCRIPT_COMMAND)
