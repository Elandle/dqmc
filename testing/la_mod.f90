module la_mod
    use iso_fortran_env, only: real32, real64, real128, output_unit, input_unit
    implicit none
    !
    ! Contains routines for working with vectors and matrices.
    ! Top performance is not guaranteed for more complicated routines
    ! (eg, no workspace queries are performed) and it is assumed that
    ! vectors are stored 1-dimensionally and matrices 2-dimensionally (eg,
    ! no leading dimension or space to skip between entries). This is the most
    ! common format for working with vectors and matrices, so it is being assumed here.
    ! The routines here can be adapted for more complex and specific situations.
    !
    ! Due to BLAS/LAPACK's use of assumed-size arrays (eg, x(*), A(lda, *)),
    ! many routines require shapes of arrays to be passed as well.
    !
    !
    integer, parameter, private :: sp = real32
    integer, parameter, private :: dp = real64
    integer, parameter, private :: ep = real128
    integer, parameter, private :: stdout = output_unit
    integer, parameter, private :: stdin = input_unit

    interface copy
        module procedure :: dvcopy
        module procedure :: dmcopy
    endinterface copy

    interface
        subroutine dcopy(n, dx, incx, dy, incy)
            import   :: dp
            integer  :: n
            real(dp) :: dx(*)
            integer  :: incx
            real(dp) :: dy(*)
            integer  :: incy
        endsubroutine dcopy
    endinterface

    contains

        subroutine dvcopy(x, y)
            !
            ! Copies:
            !
            ! x = y
            !
            ! for two vectors x and y using the BLAS routine dcopy.
            ! Assumes that x is the same length as y.
            !
            real(dp), intent(in)  :: y(:)
            real(dp), intent(out) :: x(size(y))

            call dcopy(size(y), y, 1, x, 1)


        endsubroutine dvcopy

        subroutine dmcopy(A, B, m, n, trans)
            !
            ! Copies:
            !
            ! A = B     (trans not present or trans = .false.)
            !
            ! or:
            !
            ! A = B**T  (trans = .true.)
            !
            ! for two matrices A and B using the BLAS routine dcopy.
            ! Assumes that the shape of A matches the shape of B or B**T
            ! (dependent on the value of trans)
            !
            real(dp), intent(in)           :: B(m, n)
            real(dp), intent(out)          :: A(m, n)
            integer, intent(in) :: m, n
            logical , intent(in), optional :: trans

            logical :: transpose
            integer :: j

            if (present(trans)) then
                transpose = trans
            else
                transpose = .false.
            endif

            if (.not. transpose) then
                if ((size(A, 1) .ne. size(B, 1)) .or. (size(A, 2) .ne. size(B, 2))) then
                    stop "Error in dmcopy. Copying with no transpose and matrix shapes disagree."
                endif
                call dcopy(size(B), B, 1, A, 1)
            else
                if ((size(A, 1) .ne. size(B, 2)) .or. (size(A, 2) .ne. size(B, 1))) then
                    stop "Error in dmcopy. Copying with transpose and matrix shapes disagree."
                endif
                do j = 1, size(B, 2)
                    ! A(j, :) = B(:, j)
                    call dcopy(size(B, 1), B(1, j), 1, A(j, 1), size(B, 2))
                enddo
            endif


        endsubroutine dmcopy

endmodule la_mod