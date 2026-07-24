if(NOT DEFINED EXECUTABLE OR NOT DEFINED SCENARIO OR NOT DEFINED EXPECTED_MESSAGE)
    message(FATAL_ERROR "EXECUTABLE, SCENARIO, and EXPECTED_MESSAGE are required")
endif()

execute_process(
    COMMAND "${EXECUTABLE}" --scenario "${SCENARIO}"
    RESULT_VARIABLE exit_code
    OUTPUT_VARIABLE standard_output
    ERROR_VARIABLE standard_error)

if(exit_code EQUAL 0)
    message(FATAL_ERROR "CLI unexpectedly accepted an invalid scenario")
endif()

set(combined_output "${standard_output}\n${standard_error}")
if(NOT combined_output MATCHES "${EXPECTED_MESSAGE}")
    message(FATAL_ERROR "CLI failure did not include '${EXPECTED_MESSAGE}':\n${combined_output}")
endif()
