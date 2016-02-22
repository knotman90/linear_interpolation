cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(interpolate LANGUAGES  CXX)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB SRC
    "include/*.hxx"
    "source/*.cxx"
)


add_library(interpolate ${SRC})

SET_TARGET_PROPERTIES(interpolate PROPERTIES LINKER_LANGUAGE CXX)
#target_link_libraries()
add_subdirectory(test)
