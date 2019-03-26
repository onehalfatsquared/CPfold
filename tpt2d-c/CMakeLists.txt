cmake_minimum_required(VERSION 2.8)
project(tpt2d)

add_subdirectory(support)
include_directories(support)
add_subdirectory(physics)
include_directories(physics)
add_subdirectory(exe)

file(COPY input DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/exe)
