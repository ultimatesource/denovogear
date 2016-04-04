macro(ESCAPE_STRING STR)
  string(REPLACE "\\" "\\\\" ${STR} "${${STR}}")
  string(REPLACE "\n" "\\n"  ${STR} "${${STR}}")
  string(REPLACE "\t" "\\t"  ${STR} "${${STR}}")  
endmacro()

function(CheckProcessTest PREFIX TEST)
  message(STATUS "Test ${PREFIX}.${TEST}...")

  execute_process( COMMAND ${${TEST}-CMD} 
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE  stderr
  RESULT_VARIABLE result
  WORKING_DIRECTORY "${${TEST}-WD}"
  )
  string(REPLACE "\n" "\n out> " out " out> ${stdout}")
  string(REPLACE "\n" "\n err> " err " err> ${stderr}")

  if(DEFINED ${TEST}-RESULT)
    if(NOT "${result}" STREQUAL "${${TEST}-RESULT}")
      message(FATAL_ERROR
      "Test result does not match \"${${TEST}-RESULT}\".\n"
      "Test result: ${result}\n"
      "Test output:\n"
      "${out}\n"
      "${err}" )
    endif()
  endif()

  if(DEFINED ${TEST}-RESULT-FAIL)
    if("${result}" STREQUAL "${${TEST}-RESULT-FAIL}")
      message(FATAL_ERROR
      "Test result unexpectedly matches \"${${TEST}-RESULT-FAIL}\".\n"
      "Test result: ${result}\n"
      "Test output:\n"
      "${out}\n"
      "${err}" )
    endif()
  endif()

  if(DEFINED ${TEST}-STDERR)
    foreach(test_str ${${TEST}-STDERR})
      if(NOT "${stderr}" MATCHES "${test_str}")
        ESCAPE_STRING(test_str)
        message(FATAL_ERROR
          "Test stderr does not match \"${test_str}\".\n"
          "Test result: ${result}\n"
          "Test output:\n"
          "${out}\n"
          "${err}" )
      endif()
    endforeach()
  endif()

  if(DEFINED ${TEST}-STDERR-FAIL)
    foreach(test_str ${${TEST}-STDERR-FAIL})
      if("${stderr}" MATCHES "${test_str}")
        ESCAPE_STRING(test_str)
        message(FATAL_ERROR
          "Test stderr unexpectedly matches \"${test_str}\".\n"
          "Test result: ${result}\n"
          "Test output:\n"
          "${out}\n"
          "${err}" )
      endif()
    endforeach()
  endif()

  if(DEFINED ${TEST}-STDOUT)
    foreach(test_str ${${TEST}-STDOUT})
      if(NOT "${stdout}" MATCHES "${test_str}")
        ESCAPE_STRING(test_str)
        message(FATAL_ERROR
          "Test stdout does not match \"${test_str}\".\n"
          "Test result: ${result}\n"
          "Test output:\n"
          "${out}\n"
          "${err}" )
      endif()
    endforeach()
  endif()

  if(DEFINED ${TEST}-STDOUT-FAIL)
    foreach(test_str ${${TEST}-STDOUT-FAIL})
      if("${stdout}" MATCHES "${test_str}")
        ESCAPE_STRING(test_str)
        message(FATAL_ERROR
          "Test stdout unexpectedly matches \"${test_str}\".\n"
          "Test result: ${result}\n"
          "Test output:\n"
          "${out}\n"
          "${err}" )
      endif()
    endforeach()
  endif()

endfunction()

function(CheckProcessTests PREFIX)
  foreach(TEST ${ARGN})
    CheckProcessTest("${PREFIX}" "${TEST}")
  endforeach()
endfunction()
