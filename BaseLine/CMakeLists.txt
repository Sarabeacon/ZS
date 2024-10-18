cmake_minimum_required(VERSION 3.2)

project(bclist)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

set(bclist_src ./sample.cpp ./bclist_parallel.cpp ./bclist_parallel.h ./graph.cpp ./graph.h)

add_executable(baseline ${bclist_src})


target_compile_options(baseline
    PUBLIC -I${CMAKE_SOURCE_DIR}/BaseLine/omp
    PUBLIC -fopenmp
)

target_link_options(baseline
    PUBLIC -L${CMAKE_SOURCE_DIR}/BaseLine/omp
    PUBLIC -fopenmp
)
