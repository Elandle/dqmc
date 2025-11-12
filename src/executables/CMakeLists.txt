add_executable(main
  ${CMAKE_CURRENT_SOURCE_DIR}/main.f90
)

target_link_libraries(main PRIVATE dqmc)
target_compile_options(main PRIVATE -Wall -Wextra -Wsurprising -fimplicit-none -Ofast -ffree-line-length-none -march=native)
