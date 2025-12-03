add_executable(main
  ${CMAKE_CURRENT_SOURCE_DIR}/main.f90
)

target_link_libraries(main PRIVATE dqmc)
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