cmake_minimum_required(VERSION 3.16)

project(dqmc LANGUAGES Fortran)

set(BUILD_SHARED_LIBS OFF)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options($<$<NOT:$<CONFIG:Debug>>:-flto>)
    add_link_options($<$<NOT:$<CONFIG:Debug>>:-flto>)
endif()

if(UNIX AND NOT APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()


find_package(BLAS QUIET)
find_package(LAPACK QUIET)

if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)

    find_library(
        OPENBLAS_STATIC
        NAMES openblas
        PATHS
            /usr/lib/x86_64-linux-gnu
            /usr/lib/x86_64-linux-gnu/openblas-pthread
        NO_DEFAULT_PATH
    )

    if(OPENBLAS_STATIC)

        set(BLAS_FOUND TRUE)
        set(LAPACK_FOUND TRUE)

        set(BLAS_LIBRARIES "${OPENBLAS_STATIC}")
        set(LAPACK_LIBRARIES "${OPENBLAS_STATIC}")

    else()



    endif()

endif()


if(NOT TARGET BLAS::BLAS)

    add_library(
        BLAS::BLAS
        INTERFACE
        IMPORTED
    )

    set_target_properties(
        BLAS::BLAS
        PROPERTIES
            INTERFACE_LINK_LIBRARIES
                "${BLAS_LIBRARIES};pthread"
    )
endif()

if(NOT TARGET LAPACK::LAPACK)

    add_library(
        LAPACK::LAPACK
        INTERFACE
        IMPORTED
    )

    set_target_properties(
        LAPACK::LAPACK
        PROPERTIES
            INTERFACE_LINK_LIBRARIES
                "${LAPACK_LIBRARIES}"
    )

endif()

add_subdirectory(src)