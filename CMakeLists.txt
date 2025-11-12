cmake_minimum_required(VERSION 3.20)

project(dqmc LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_subdirectory(src)