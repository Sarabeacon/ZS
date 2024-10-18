cmake_minimum_required(VERSION 3.16)

project(pqclique)

add_subdirectory( ${CMAKE_SOURCE_DIR}/BaseLine )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
include_directories(include)

file(GLOB demo_SORCE_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")
add_executable(ZS ${demo_SORCE_FILES})

target_link_options(ZS 
    PUBLIC -lpthread
)



