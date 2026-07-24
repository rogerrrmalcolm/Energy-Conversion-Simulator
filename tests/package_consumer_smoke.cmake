foreach(required_variable PROJECT_BINARY_DIR PACKAGE_SOURCE_DIR TEST_BINARY_DIR EXECUTABLE_SUFFIX)
    if(NOT DEFINED ${required_variable})
        message(FATAL_ERROR "${required_variable} is required")
    endif()
endforeach()

function(run_checked)
    execute_process(
        COMMAND ${ARGV}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error)
    if(NOT ${result} EQUAL 0)
        message(FATAL_ERROR "Command failed: ${ARGV}\n${output}\n${error}")
    endif()
endfunction()

set(install_prefix "${TEST_BINARY_DIR}/install")
set(consumer_build "${TEST_BINARY_DIR}/consumer")
file(REMOVE_RECURSE "${TEST_BINARY_DIR}")
file(MAKE_DIRECTORY "${TEST_BINARY_DIR}")

set(install_command "${CMAKE_COMMAND}" --install "${PROJECT_BINARY_DIR}" --prefix "${install_prefix}")
if(DEFINED CONFIG AND NOT CONFIG STREQUAL "")
    list(APPEND install_command --config "${CONFIG}")
endif()
run_checked(${install_command})

set(configure_command
    "${CMAKE_COMMAND}"
    -S "${PACKAGE_SOURCE_DIR}"
    -B "${consumer_build}"
    "-DCMAKE_PREFIX_PATH=${install_prefix}")
if(DEFINED CONFIG AND NOT CONFIG STREQUAL "")
    list(APPEND configure_command "-DCMAKE_BUILD_TYPE=${CONFIG}")
endif()
run_checked(${configure_command})

set(build_command "${CMAKE_COMMAND}" --build "${consumer_build}")
if(DEFINED CONFIG AND NOT CONFIG STREQUAL "")
    list(APPEND build_command --config "${CONFIG}")
endif()
run_checked(${build_command})

if(DEFINED CONFIG AND NOT CONFIG STREQUAL "")
    set(configured_consumer_executable
        "${consumer_build}/${CONFIG}/package_consumer${EXECUTABLE_SUFFIX}")
else()
    set(configured_consumer_executable "")
endif()
set(consumer_executable "${consumer_build}/package_consumer${EXECUTABLE_SUFFIX}")
if(NOT configured_consumer_executable STREQUAL "" AND
   EXISTS "${configured_consumer_executable}")
    set(consumer_executable "${configured_consumer_executable}")
endif()
run_checked("${consumer_executable}")
