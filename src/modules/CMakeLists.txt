
add_library(dqmc
      ${CMAKE_CURRENT_SOURCE_DIR}/blas_interface_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/bmult_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/checkerboard_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/convenientla_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/customla_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/equalgreens_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/expm_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/lapack_interface_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/measurements_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/metropolisratios_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/numbertypes.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/polarization_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/printing_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/readinputfile_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/simulate_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/simulationsetup_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/statistics_mod.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/utilities_mod.f90
)

set(module_directory "${CMAKE_CURRENT_BINARY_DIR}/module_directory")

set_target_properties(dqmc PROPERTIES
      Fortran_MODULE_DIRECTORY "${module_directory}"
      OUTPUT_NAME dqmc
)

target_include_directories(dqmc
  PUBLIC
    "${module_directory}"
)

target_link_libraries(dqmc
    PUBLIC
        LAPACK::LAPACK
        BLAS::BLAS
)

target_compile_options(dqmc PRIVATE -Wall -Wextra -Wsurprising -fimplicit-none -Ofast -ffree-line-length-none -march=native)

install(TARGETS dqmc
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
