cmake_minimum_required(VERSION 3.20)
project(black_hole_energy_simulation VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BH_BUILD_TESTS "Build the black-hole model tests" ON)
option(BH_ENABLE_AVX2 "Build the four-lane AVX2 numerical kernels" OFF)

function(bh_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endfunction()

add_library(black_hole_models
    src/physics/algebraic_model.cpp
    src/physics/penrose_model.cpp
    src/io/penrose_scenario_io.cpp
    src/optimization/dijkstra.cpp
    src/integrators/schwarzschild_geodesic.cpp
    src/integrators/kerr_geodesic.cpp
    src/physics/plasma_model.cpp)
add_library(bh::models ALIAS black_hole_models)
target_compile_features(black_hole_models PUBLIC cxx_std_20)
set_target_properties(black_hole_models PROPERTIES
    EXPORT_NAME models)
target_include_directories(black_hole_models
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(BH_ENABLE_AVX2)
    target_compile_definitions(black_hole_models PRIVATE BH_ENABLE_AVX2=1)
    if(MSVC)
        target_compile_options(black_hole_models PRIVATE /arch:AVX2)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(black_hole_models PRIVATE -mavx2)
    else()
        message(FATAL_ERROR "BH_ENABLE_AVX2 requires an MSVC, GCC, or Clang x86 compiler")
    endif()
endif()

bh_enable_warnings(black_hole_models)

add_executable(black_hole_demo
    src/app/main.cpp
    src/app/cli_json.cpp)
target_link_libraries(black_hole_demo PRIVATE bh::models)
target_compile_definitions(black_hole_demo PRIVATE BH_VERSION="${PROJECT_VERSION}")
bh_enable_warnings(black_hole_demo)
set_target_properties(black_hole_demo PROPERTIES OUTPUT_NAME black-hole-sim)

if(BH_BUILD_TESTS)
    enable_testing()
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_executable(model_tests tests/model_tests.cpp)
    target_link_libraries(model_tests PRIVATE bh::models)
    target_compile_definitions(model_tests PRIVATE BH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    bh_enable_warnings(model_tests)
    add_test(NAME model_tests COMMAND model_tests)

    add_test(NAME cli_help COMMAND black_hole_demo --help)
    set_tests_properties(cli_help PROPERTIES
        PASS_REGULAR_EXPRESSION "Usage:")

    add_test(NAME cli_version COMMAND black_hole_demo --version)
    set_tests_properties(cli_version PROPERTIES
        PASS_REGULAR_EXPRESSION "black-hole-sim ${PROJECT_VERSION}")

    set(BH_INTERACTIVE_TEST_ARGUMENTS
        --executable $<TARGET_FILE:black_hole_demo>
        --input "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/interactive_pipeline.txt")
    if(BH_ENABLE_AVX2)
        list(APPEND BH_INTERACTIVE_TEST_ARGUMENTS --require-avx2)
    endif()
    add_test(NAME cli_sequential_interactive COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/assert_interactive_pipeline.py"
        ${BH_INTERACTIVE_TEST_ARGUMENTS})

    add_test(NAME cli_algebraic COMMAND black_hole_demo algebraic 1.98847e31 0.9)
    set_tests_properties(cli_algebraic PROPERTIES
        PASS_REGULAR_EXPRESSION "rotational reservoir:")

    add_test(NAME cli_algebraic_range COMMAND black_hole_demo algebraic-range
        1.590776e31 1.98847e31 2.386164e31 0.8 0.9 0.99)
    set_tests_properties(cli_algebraic_range PROPERTIES
        PASS_REGULAR_EXPRESSION "upper uncertainty:")

    add_test(NAME cli_toy_plasma COMMAND black_hole_demo toy-plasma 1.0 1e-8 10.0 0.9 2.0)
    set_tests_properties(cli_toy_plasma PROPERTIES
        PASS_REGULAR_EXPRESSION "Reduced toy-plasma transport estimate")

    add_test(NAME cli_penrose_reference COMMAND black_hole_demo scenario
        "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/equatorial_penrose_reference.cfg")
    set_tests_properties(cli_penrose_reference PROPERTIES
        PASS_REGULAR_EXPRESSION "status: physically_feasible"
        FAIL_REGULAR_EXPRESSION "Integration diagnostics")

    add_test(NAME cli_penrose_reference_summary COMMAND black_hole_demo scenario
        "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/equatorial_penrose_reference.cfg" --summary)
    set_tests_properties(cli_penrose_reference_summary PROPERTIES
        PASS_REGULAR_EXPRESSION "status: physically_feasible"
        FAIL_REGULAR_EXPRESSION "Integration diagnostics")

    add_test(NAME cli_penrose_reference_verbose COMMAND black_hole_demo scenario
        "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/equatorial_penrose_reference.cfg" --verbose)
    set_tests_properties(cli_penrose_reference_verbose PROPERTIES
        PASS_REGULAR_EXPRESSION "Integration diagnostics")

    add_test(NAME cli_penrose_dijkstra_goal_fixture COMMAND black_hole_demo search
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/equatorial_penrose_dijkstra_goal.cfg")
    set_tests_properties(cli_penrose_dijkstra_goal_fixture PROPERTIES
        PASS_REGULAR_EXPRESSION "found_goal")

    add_test(NAME cli_penrose_dijkstra_verbose COMMAND black_hole_demo search
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/equatorial_penrose_dijkstra_goal.cfg"
        --verbose --no-progress)
    set_tests_properties(cli_penrose_dijkstra_verbose PROPERTIES
        PASS_REGULAR_EXPRESSION "Integration diagnostics")

    add_test(NAME cli_penrose_phase_map COMMAND black_hole_demo map
        "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/equatorial_penrose_dijkstra_smoke.cfg")
    set_tests_properties(cli_penrose_phase_map PROPERTIES
        PASS_REGULAR_EXPRESSION "completed /.*single-thread")

    add_test(NAME cli_penrose_phase_map_verbose COMMAND black_hole_demo map
        "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/equatorial_penrose_dijkstra_smoke.cfg"
        --verbose --no-progress)
    set_tests_properties(cli_penrose_phase_map_verbose PROPERTIES
        PASS_REGULAR_EXPRESSION "Integration diagnostics")

    if(BH_ENABLE_AVX2)
        add_test(NAME cli_penrose_avx2_kernel COMMAND black_hole_demo map
            "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/equatorial_penrose_avx2_smoke.cfg")
        set_tests_properties(cli_penrose_avx2_kernel PROPERTIES
            PASS_REGULAR_EXPRESSION "AVX2 batches: +[1-9][0-9]*")
    endif()

    add_test(NAME cli_penrose_complete_window_fallback COMMAND black_hole_demo search
        "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/equatorial_penrose_dijkstra_smoke.cfg")
    set_tests_properties(cli_penrose_complete_window_fallback PROPERTIES
        PASS_REGULAR_EXPRESSION "target not reached; showing the best validated candidate")

    add_test(NAME cli_legacy_flag_alias COMMAND black_hole_demo --algebraic 1.98847e31 0.9)
    set_tests_properties(cli_legacy_flag_alias PROPERTIES
        PASS_REGULAR_EXPRESSION "rotational reservoir:")

    add_test(NAME cli_json_algebraic COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/assert_cli_json.py"
        --executable $<TARGET_FILE:black_hole_demo>
        --expected-command algebraic -- algebraic 1.98847e31 0.9)

    add_test(NAME cli_json_search COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/assert_cli_json.py"
        --executable $<TARGET_FILE:black_hole_demo>
        --expected-command search --expected-status found_goal -- search
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/equatorial_penrose_dijkstra_goal.cfg")

    if(BH_ENABLE_AVX2)
        add_test(NAME cli_json_avx2_map COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/assert_cli_json.py"
            --executable $<TARGET_FILE:black_hole_demo>
            --expected-command map --expected-status completed --require-avx2
            --require-progress -- map
            "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/equatorial_penrose_avx2_smoke.cfg")
    endif()

    add_test(NAME cli_invalid_penrose_scenario COMMAND ${CMAKE_COMMAND}
        "-DEXECUTABLE=$<TARGET_FILE:black_hole_demo>"
        "-DSCENARIO=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/invalid_penrose_unknown_key.cfg"
        "-DEXPECTED_MESSAGE=unknown scenario key"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/assert_cli_failure.cmake")

    add_test(NAME installed_package_consumer COMMAND ${CMAKE_COMMAND}
        "-DPROJECT_BINARY_DIR=${CMAKE_BINARY_DIR}"
        "-DPACKAGE_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests/package_consumer"
        "-DTEST_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/package-consumer-smoke"
        "-DCONFIG=$<CONFIG>"
        "-DEXECUTABLE_SUFFIX=${CMAKE_EXECUTABLE_SUFFIX}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/package_consumer_smoke.cmake")
endif()

install(TARGETS black_hole_models
    EXPORT BlackHoleModelsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS black_hole_demo
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY scenarios/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/black_hole_energy_simulation/scenarios)
install(EXPORT BlackHoleModelsTargets
    FILE BlackHoleModelsTargets.cmake
    NAMESPACE bh::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BlackHoleModels)

configure_package_config_file(
    cmake/BlackHoleModelsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/BlackHoleModelsConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BlackHoleModels)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/BlackHoleModelsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/BlackHoleModelsConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/BlackHoleModelsConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BlackHoleModels)
