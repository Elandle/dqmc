module qr_mod
    !
    ! Module for calling QR factorisations conveniently.
    ! Focus is on convenience of calling for testing things out, not performance.
    !
    ! At the moment these only work for square (m = n) matrices.
    ! Todo: fix for m =/= n matrices
    !
    use iso_fortran_env, only: real64, output_unit
    implicit none

    integer, parameter, private :: dp       = real64
    integer, parameter, private :: terminal = output_unit

    interface
        subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
            INTEGER            INFO, LDA, LWORK, M, N
            DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        endsubroutine dgeqrf

        subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
            INTEGER            INFO, K, LDA, LWORK, M, N
            DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        endsubroutine dorgqr

        subroutine dgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info)
            INTEGER            INFO, LDA, LWORK, M, N
            INTEGER            JPVT( * )
            DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        endsubroutine dgeqp3
    end interface

    contains

        subroutine qrp(A, P, Q, R)
            !
            ! A: m x n
            ! P: n x n
            ! Q: m x m
            ! R: m x n
            !
            ! AP = QR
            !
            ! Input arguments
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: P(size(A, 2), size(A, 2))
            real(dp), intent(out) :: Q(size(A, 1), size(A, 1))
            real(dp), intent(out) :: R(size(A, 1), size(A, 2))

            ! Subroutine variables
            real(dp) :: matrix(size(A, 1), size(A, 2))
            integer  :: m, n

            ! Extra variables for LAPACK
            real(dp) :: tau(min(size(A, 1), size(A, 2)))
            real(dp) :: work(10*size(A, 2))
            integer  :: jpvt(size(A, 2))
            integer  :: lwork
            integer  :: info

            m = size(A, 1); n = size(A, 2)
            lwork = size(work)
            info = 0
            tau = 0.0_dp
            work = 0.0_dp
            jpvt = 0

            matrix = A


            call dgeqp3(m, n, matrix, m, jpvt, tau, work, lwork, info)
            if (info .lt. 0) then
                write(terminal, "(a, i4)") "Error in subroutine qrp calling dgeqp3. LAPACK info = ", info
            endif
            P = jpvt_to_p(jpvt)
            Q = qr_to_q(matrix, tau)
            R = triu(matrix)
        endsubroutine qrp

        subroutine qr(A, Q, R)
            !
            ! A: m x n
            ! Q: m x m
            ! R: m x n
            !
            ! A = QR
            !
            ! Input arguments
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: Q(size(A, 1), size(A, 1))
            real(dp), intent(out) :: R(size(A, 1), size(A, 2))

            ! Subroutine variables
            real(dp) :: matrix(size(A, 1), size(A, 2))
            integer  :: m, n

            ! Extra variables for LAPACK
            real(dp) :: tau(min(size(A, 1), size(A, 2)))
            real(dp) :: work(10*size(A, 2))
            integer  :: lwork
            integer  :: info


            m = size(A, 1); n = size(A, 2)
            lwork = size(work)
            info = 0
            tau = 0.0_dp
            work = 0.0_dp

            matrix = A

            call dgeqrf(m, n, matrix, m, tau, work, lwork, info)
            if (info .lt. 0) then
                write(terminal, "(a, i4)") "Error in subroutine qr calling dgeqrf. LAPACK info = ", info
            endif
            Q = qr_to_q(matrix, tau)
            R = triu(matrix)
        endsubroutine qr

        function triu(A, k) result(U)
            !
            ! Returns a matrix the same shape as A with the elements below the kth diagonal zeroed.
            ! If not present, k = 0 by default. Then the upper triangular part of A is returned.
            ! k can be positive, negative, or zero.
            !
            real(dp),           intent(in) :: A(:, :)
            integer , optional, intent(in) :: k
            real(dp)                       :: U(size(A, 1), size(A, 2))

            integer :: i, j, k0

            if (present(k)) then
                k0 = k
            else
                k0 = 0
            endif

            U = 0.0_dp
            do i = 1, size(A, 1)
                do j = max(1, i+k0), size(A, 2)
                    U(i, j) = A(i, j)
                enddo
            enddo
        endfunction triu

        function qr_to_q(A, tau) result(Q)
            real(dp), intent(in) :: A(:, :)
            real(dp), intent(in) :: tau(min(size(A, 1), size(A, 2)))

            real(dp) :: Q(size(A, 1), size(A, 1))
            real(dp) :: matrix(size(A, 1), size(A, 2))
            real(dp) :: work(size(A, 2) * 256)
            integer  :: m, n, k, lwork, info

            matrix = A
            m = size(A, 1)
            n = size(A, 2)
            k = min(m, n)
            lwork = size(work)

            call dorgqr(m, n, k, matrix, m, tau, work, lwork, info)

            Q = matrix(:, :m)
        endfunction qr_to_q

        function jpvt_to_p(jpvt) result(P)
            integer , intent(in) :: jpvt(:)
            real(dp)             :: P(size(jpvt), size(jpvt))

            integer :: i

            P = 0.0_dp
            do i = 1, size(jpvt)
                P(jpvt(i), i) = 1.0_dp
            enddo
        endfunction jpvt_to_p
endmodule qr_mod