module numbertypes
    use iso_fortran_env, only: real32, real64
    implicit none
    !
    ! Specifies number types for portability
    !
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64

endmodule numbertypes