cmake_minimum_required(VERSION 3.20)
project(black_hole_models_package_smoke_test LANGUAGES CXX)

find_package(BlackHoleModels 0.1 CONFIG REQUIRED)

add_executable(package_consumer main.cpp)
target_link_libraries(package_consumer PRIVATE bh::models)
