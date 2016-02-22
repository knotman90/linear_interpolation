cmake_minimum_required(VERSION 3.1 FATAL_ERROR)


find_package(GTest REQUIRED)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)


include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SRC
    "*.hxx"
    "*.cxx"
)

add_executable(tests ${SRC})
#target_compile_features(tests PRIVATE cxx_deleted_functions)
set_property(TARGET tests PROPERTY CXX_STANDARD 14)
target_link_libraries(tests ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} pthread)
