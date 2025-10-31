program main
    use iso_fortran_env, only: dp => real64, terminal => output_unit
    implicit none

    integer, parameter :: m = 2
    
    real(dp) :: A(m, m)

    A = id(m)
    call print_matrix(A, terminal, "A = ")
    call sub1(A)
    call print_matrix(A, terminal, "A = ")

    



    contains

        subroutine sub1(A)
            real(dp), intent(out) :: A(:, :)

            A = 2.0_dp * id(size(A, 1))

        endsubroutine sub1

        subroutine sub2(A)
            real(dp), intent(out) :: A(:, :)

            A = 5.0_dp * id(size(A, 1))

        endsubroutine sub2


        function id(m) result(A)
            integer, intent(in) :: m
            
            real(dp) :: A(m, m)
            integer  :: i

            A = 0.0_dp
            do i = 1, m
                A(i, i) = 1.0_dp
            enddo


        endfunction id


        subroutine print_matrix(A, ounit, message)
            real(dp)        , intent(in)           :: A(:, :)
            integer         , intent(in)           :: ounit
            character(len=*), intent(in), optional :: message
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            if (present(message)) then
                write (ounit, "(a)") message
            endif

            do i = 1, m
                do j = 1, n
                    write(ounit, "(f17.8)", advance="no") A(i, j)
                enddo
                write(ounit, "(a)") ""
            enddo


        endsubroutine print_matrix




endprogram main