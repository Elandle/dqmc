
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
      ${CMAKE_CURRENT_SOURCE_DIR}/stduse.f90
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

# -march=native
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(dqmc PRIVATE
        -Wall -Wextra -Wsurprising       
        -fimplicit-none                  
        -ffree-line-length-none          
        -Ofast
        -static                                             
    )
# -xHost
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    target_compile_options(dqmc PRIVATE
        -warn all                                          
        -extend_source                   
        -Ofast
        -static-intel
        -static                                                     
    )
endif()


install(TARGETS dqmc
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
