program main
    use iso_fortran_env, only: real64, output_unit
    use la_mod
    use printing_mod
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: stdout = output_unit
    integer, parameter :: m = 4, n = 4

    real(dp) :: A(m, n), B(m, n)

    !call random_number(B)

    !call print_matrix(B, stdout, "B = ")

    !call copy(A, b)

    !call print_matrix(A, stdout, "A = ")





endprogram main