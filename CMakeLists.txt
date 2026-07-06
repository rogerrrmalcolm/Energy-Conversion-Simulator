cmake_minimum_required(VERSION 3.20)
project(black_hole_energy_simulation VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(BH_BUILD_TESTS "Build the black-hole model tests" ON)

add_library(black_hole_models
    src/physics/algebraic_model.cpp
    src/integrators/schwarzschild_geodesic.cpp
    src/integrators/kerr_geodesic.cpp
    src/physics/plasma_model.cpp)
add_library(bh::models ALIAS black_hole_models)
target_compile_features(black_hole_models PUBLIC cxx_std_20)
set_target_properties(black_hole_models PROPERTIES
    CXX_EXTENSIONS OFF
    EXPORT_NAME models)
target_include_directories(black_hole_models
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(MSVC)
    target_compile_options(black_hole_models PRIVATE /W4 /permissive-)
else()
    target_compile_options(black_hole_models PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(black_hole_demo src/app/main.cpp)
target_link_libraries(black_hole_demo PRIVATE bh::models)

if(BH_BUILD_TESTS)
    enable_testing()
    add_executable(model_tests tests/model_tests.cpp)
    target_link_libraries(model_tests PRIVATE bh::models)
    add_test(NAME model_tests COMMAND model_tests)
endif()

install(TARGETS black_hole_models
    EXPORT BlackHoleModelsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
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
