module convenientla_mod
    use numbertypes
    use customla_mod
    implicit none

    ! Contains convenient, not high-performing, linear algebra procedures.
    ! Mainly for testing
    !
    ! For example, the update:
    !
    !     A = A * B
    !
    ! without having to supply a workspace
    !

    contains


        subroutine right_matrixmultiply(A, B, n)
            !
            ! Updates:
            !
            ! A = A * B
            !
            ! where A and B are n x n matrices
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: B(n, n)
            integer , intent(in)    :: n

            real(dp), allocatable :: C(:, :)

            allocate(C(n, n))

            ! C = A
            call copy_matrix(A, C, n)
            ! A = C * B
            call dgemm('n', 'n', n, n, n, 1.0_dp, C, n, B, n, 0.0_dp, A, n)


        endsubroutine


        subroutine left_matrixmultiply(A, B, n)
            !
            ! Updates:
            !
            ! A = B * A
            !
            ! where A and B are n x n matrices
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: B(n, n)
            integer , intent(in)    :: n

            real(dp), allocatable :: C(:, :)

            allocate(C(n, n))

            ! C = A
            call copy_matrix(A, C, n)
            ! A = B * C
            call dgemm('n', 'n', n, n, n, 1.0_dp, B, n, C, n, 0.0_dp, A, n)


        endsubroutine


        real(dp) function twonorm(A) result(norm)
            !
            ! Returns the 2 norm of the matrix A
            !
            real(dp), intent(in)  :: A(:, :)
            real(dp), allocatable :: B(:, :)
            real(dp), allocatable :: S(:)
            real(dp), allocatable :: work(:)
            integer               :: m
            integer               :: n
            integer               :: lwork
            integer               :: info
            
            m = size(A, 1)
            n = size(A, 2)
            lwork = 5 * (3 * min(m, n) + max(m, n) + 5 * min(m, n))

            allocate(B(m, n))
            allocate(S(min(m, n)))
            allocate(work(lwork))

            call dlacpy('A', m, n, A, m, B, m)
            call dgesvd('N', 'N', m, n, B, m, S, work, lwork, work, lwork, work, lwork, info)

            norm = S(1)

            
        endfunction twonorm






endmodule convenientla_mod