cmake_minimum_required(VERSION 3.16)

project(dqmc LANGUAGES Fortran)

set(BUILD_SHARED_LIBS OFF)

if(UNIX AND NOT APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    set_target_properties(
        BLAS::BLAS
        PROPERTIES
            INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES};pthread"
    )
endif()

if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_target_properties(
        LAPACK::LAPACK
        PROPERTIES
            INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}"
    )
endif()

add_subdirectory(src)
