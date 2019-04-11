cmake_minimum_required(VERSION 2.8)
project(tpt2d)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp ")

add_subdirectory(nauty)
include_directories(nauty)
add_subdirectory(support)
include_directories(support)
add_subdirectory(physics)
include_directories(physics)
add_subdirectory(tpt)
include_directories(tpt)
add_subdirectory(exe)

file(COPY input DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/exe)
