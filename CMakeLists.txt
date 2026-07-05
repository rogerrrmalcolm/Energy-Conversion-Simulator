cmake_minimum_required(VERSION 3.20)
project(black_hole_energy_simulation VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(black_hole_models
    src/physics/algebraic_model.cpp
    src/integrators/schwarzschild_geodesic.cpp
    src/integrators/kerr_geodesic.cpp
    src/physics/plasma_model.cpp)
target_include_directories(black_hole_models PUBLIC include)

if(MSVC)
    target_compile_options(black_hole_models PRIVATE /W4 /permissive-)
else()
    target_compile_options(black_hole_models PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(black_hole_demo src/app/main.cpp)
target_link_libraries(black_hole_demo PRIVATE black_hole_models)

enable_testing()
add_executable(model_tests tests/model_tests.cpp)
target_link_libraries(model_tests PRIVATE black_hole_models)
add_test(NAME model_tests COMMAND model_tests)
