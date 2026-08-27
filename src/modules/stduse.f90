module stduse
    use iso_fortran_env, only: real32, real64, input_unit, output_unit, iostat_end
    implicit none

    integer, parameter :: sp     = real32
    integer, parameter :: dp     = real64
    integer, parameter :: stdin  = input_unit
    integer, parameter :: stdout = output_unit
endmodule stduse