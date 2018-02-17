macro(ESCAPE_STRING STR)
  string(REPLACE ";" "%3B" ${STR} "${${STR}}")
  string(REPLACE "\\." "%2E" ${STR} "${${STR}}")
  string(REPLACE "\\" "\\\\" ${STR} "${${STR}}")
  string(REPLACE "\n" "\\n"  ${STR} "${${STR}}")
  string(REPLACE "\r" "\\r"  ${STR} "${${STR}}")
  string(REPLACE "\t" "\\t"  ${STR} "${${STR}}")  
endmacro()

set(failures)
set(test_failed FALSE)

function(ReportFailure MSG)
  set(failures ${failures} ${MSG} PARENT_SCOPE)
  set(test_failed TRUE PARENT_SCOPE)
endfunction()

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
      ReportFailure("Test result does not match \"${${TEST}-RESULT}\".")
    endif()
  endif()

  if(DEFINED ${TEST}-RESULT-FAIL)
    if("${result}" STREQUAL "${${TEST}-RESULT-FAIL}")
      ReportFailure("Test result unexpectedly matches \"${${TEST}-RESULT-FAIL}\".")
    endif()
  endif()

  if(DEFINED ${TEST}-STDERR)
    foreach(test_str ${${TEST}-STDERR})
      if(NOT "${stderr}" MATCHES "${test_str}")
        ESCAPE_STRING(test_str)
        ReportFailure("Test stderr does not match \"${test_str}\".")
      endif()
    endforeach()
  endif()

  if(DEFINED ${TEST}-STDERR-FAIL)
    foreach(test_str ${${TEST}-STDERR-FAIL})
      if("${stderr}" MATCHES "${test_str}")
        ESCAPE_STRING(test_str)
        ReportFailure("Test stderr unexpectedly matches \"${test_str}\".")
      endif()
    endforeach()
  endif()

  if(DEFINED ${TEST}-STDOUT)
    foreach(test_str ${${TEST}-STDOUT})
      if(NOT "${stdout}" MATCHES "${test_str}")
        ESCAPE_STRING(test_str)
        ReportFailure("Test stdout does not match \"${test_str}\".")
      endif()
    endforeach()
  endif()

  if(DEFINED ${TEST}-STDOUT-FAIL)
    foreach(test_str ${${TEST}-STDOUT-FAIL})
      if("${stdout}" MATCHES "${test_str}")
        ESCAPE_STRING(test_str)
        ReportFailure("Test stdout unexpectedly matches \"${test_str}\".")
      endif()
    endforeach()
  endif()

  if(test_failed)
    message(
      " res> ${result}\n"
      "${out}\n"
      "${err}"
    )
    foreach(msg ${failures})
      string(REPLACE "%3B" ";" msg "${msg}")
      string(REPLACE "%2E" "." msg "${msg}")
      message("  FAILURE: ${msg}")
    endforeach()
    message(SEND_ERROR "Test failed.")
  endif()
endfunction()

function(CheckProcessTests PREFIX)
  foreach(TEST ${ARGN})
    CheckProcessTest("${PREFIX}" "${TEST}")
  endforeach()
endfunction()
