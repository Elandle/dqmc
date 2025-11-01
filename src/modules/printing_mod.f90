module printing_mod
    use iso_fortran_env, only: real64
    implicit none

    integer, parameter, private :: dp = real64

    public :: print_matrix
    public :: print_vector

    contains

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
                    ! write(ounit, "(e17.8)", advance="no") A(i, j)
                enddo
                write(ounit, "(a)") ""
            enddo
        endsubroutine print_matrix

        subroutine print_vector(v, ounit, advance)
            real(dp), intent(in) :: v(:)
            integer , intent(in) :: ounit
            logical , optional   :: advance
            
            integer :: m, n, i, j
            logical :: adv

            if (present(advance)) then
                adv = advance
            else
                adv = .false.
            endif

            m = size(v)

            if (adv) then
                do i = 1, m
                    write(ounit, "(f17.8)") v(i)
                enddo
            else
                do i = 1, m
                    write(ounit, "(f17.8)", advance="no") v(i)
                enddo
                write(ounit, "(a)") ""
            endif
        endsubroutine print_vector

        subroutine print_integer_matrix(A, ounit, message)
            integer         , intent(in)           :: A(:, :)
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
                    write(ounit, "(i6)", advance="no") A(i, j)
                    ! write(ounit, "(e17.8)", advance="no") A(i, j)
                enddo
                write(ounit, "(a)") ""
            enddo
        endsubroutine print_integer_matrix
        
endmodule printing_mod