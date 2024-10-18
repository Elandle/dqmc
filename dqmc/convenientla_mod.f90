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


        subroutine exponentiate(A)
            !
            ! Updates:
            !
            ! A = exp(A)
            !
            ! Implementation is coded just to work: performance is not good at all
            ! To avoid working with LAPACK's real and imaginary components in
            ! diagonalizing a real matrix, A is first copied to a complex matrix,
            ! which is then diagonalized, exponentiated, and copied back to A
            !
            real(dp), intent(inout) :: A(:, :)

            integer :: m, n, info, i, j
            complex(dp), allocatable :: Z(:, :), W(:), VL(:, :), VR(:, :), work(:), VRinv(:, :)
            real(dp) , allocatable :: rwork(:)
            integer, allocatable :: P(:)

            m = size(A, 1)
            n = size(A, 2)

            allocate(Z(n, n), W(n), VL(n, n), VR(n, n), work(10*n), rwork(2*n), VRinv(n, n))
            allocate(P(n))

            call zlacp2('a', n, n, A, n, Z, n)
            call zgeev('V', 'V', n, Z, n, W, VL, n, VR, n, work, 10*n, rwork, info)
            
            W = exp(W)

            call zlacpy('a', n, n, VR, n, VRinv, n)
            P = 0
            call zgetrf(n, n, VRinv, n, P, info)
            call zgetri(n, VRinv, n, P, work, 10*n, info)


            do i = 1, n
                call zscal(n, W(i), VR(1, i), 1)
            enddo

            call zgemm('n', 'n', n, n, n, cmplx(1.0_dp, kind=dp), VR, n, VRinv, n, cmplx(0.0_dp, kind=dp), Z, n)

            A = real(Z, kind=dp)





        endsubroutine exponentiate









endmodule convenientla_mod