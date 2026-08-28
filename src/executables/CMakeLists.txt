add_executable(
    main
    ${CMAKE_CURRENT_SOURCE_DIR}/main.f90
)

target_link_libraries(
    main
    PRIVATE
        dqmc
)

set(
    optimization
    "nativehigh"
    CACHE
    STRING
    "Optimization level"
)

set_property(
    CACHE optimization
    PROPERTY STRINGS
        nativehigh
        puma
        portable
        debug
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

    target_link_options(
        main
        PRIVATE
            -static
    )

    target_compile_options(
        dqmc
        PRIVATE
            -Wall
            -Wextra
            -Wsurprising
            -fimplicit-none
            -ffree-line-length-none
    )

    if(optimization STREQUAL "nativehigh")
        target_compile_options(
            dqmc
            PRIVATE
                -Ofast
                -march=native
                -mtune=native
        )

    elseif(optimization STREQUAL "puma")

        target_compile_options(
            dqmc
            PRIVATE
                -Ofast
                -march=sandybridge
                -mtune=generic
        )

    elseif(optimization STREQUAL "portable")

        target_compile_options(
            dqmc
            PRIVATE
                -O3
                -march=x86-64-v2
                -mtune=generic
        )

    elseif(optimization STREQUAL "debug")

        target_compile_options(
            dqmc
            PRIVATE
                -O0
                -g
                -fcheck=all
                -fbacktrace
        )

    else()

endif()