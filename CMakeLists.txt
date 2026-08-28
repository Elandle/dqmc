cmake_minimum_required(VERSION 3.16)

project(dqmc LANGUAGES Fortran)

set(BUILD_SHARED_LIBS OFF)

if(UNIX AND NOT APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_subdirectory(src)