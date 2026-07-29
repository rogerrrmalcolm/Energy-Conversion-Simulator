if(NOT DEFINED EXECUTABLE OR NOT DEFINED INPUT_FILE OR NOT DEFINED REQUIRED_MESSAGES)
    message(FATAL_ERROR "EXECUTABLE, INPUT_FILE, and REQUIRED_MESSAGES are required")
endif()

execute_process(
    COMMAND "${EXECUTABLE}" --interactive
    INPUT_FILE "${INPUT_FILE}"
    RESULT_VARIABLE exit_code
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error_output)

if(NOT exit_code EQUAL 0)
    message(FATAL_ERROR
        "Interactive CLI failed with exit code ${exit_code}.\nstdout:\n${output}\nstderr:\n${error_output}")
endif()

string(REPLACE "|" ";" required_message_list "${REQUIRED_MESSAGES}")
foreach(required_message IN LISTS required_message_list)
    string(FIND "${output}" "${required_message}" expected_position)
    if(expected_position EQUAL -1)
        message(FATAL_ERROR
            "Interactive CLI output did not contain '${required_message}'.\nstdout:\n${output}")
    endif()
endforeach()
