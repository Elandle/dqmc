program main
    use iso_fortran_env, only: dp => real64, terminal => output_unit
    implicit none


    integer , parameter :: m = 8, n = 8
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)


    real(dp)              :: A(m, n), P(n, n), Q(n, n)
    real(dp), allocatable :: R(:, :)
    integer  :: k, pp

    call random_number(A)

    k  = 4
    pp = 0
    allocate(R(k, n))
    call ssrqrcp(A, P)

    call print_matrix(P, terminal, "P = ")


    contains


        subroutine ssrqrcp(A, P)
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: P(size(A, 2), size(A, 2))

            call qrcp(A, P)
            call print_matrix(P, terminal, "P in ssrqrcp")

        endsubroutine ssrqrcp


        subroutine qrcp(A, P)
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: P(size(A, 2), size(A, 2))


            integer  :: jpvt(size(A, 2))

            jpvt = 0
            P = 0.0_dp

            jpvt(2) = 3
            P = jpvt_to_p(jpvt)

            
        endsubroutine qrcp


        function dgeq_to_r(A) result(R)
            real(dp), intent(in) :: A(:, :)

            real(dp) :: R(size(A, 1), size(A, 2))
            integer  :: i, j

            R = 0.0_dp
            do i = 1, size(R, 1)
                do j = i, size(R, 2)
                    R(i, j) = A(i, j)
                enddo
            enddo



        endfunction dgeq_to_r



        function dgeq_to_q2(A, tau) result(Q)
            real(dp), intent(in) :: A(:, :)
            real(dp), intent(in) :: tau(min(size(A, 1), size(A, 2)))

            real(dp) :: Q(min(size(A, 1), size(A, 2)), min(size(A, 1), size(A, 2)))
            real(dp) :: matrix(size(A, 1), size(A, 2))
            integer  :: m, n, k
            real(dp) :: work(128*size(A, 1))
            integer  :: lwork
            integer  :: info

            m = size(A, 1); n = size(A, 2)
            lwork = size(work)
            info = 0
            work = 0.0_dp
            k = min(m, n)
            matrix = A
            Q = 0.0_dp

            print *, "m = ", m, "n = ", n
            print *, "matrix dimensions = ", size(matrix, 1), size(matrix, 2)

            call dorgqr(k, k, k, matrix, m, tau, work, lwork, info)
            if (info .lt. 0) then
                write(terminal, "(a, i4)") "Error in subroutine dgeq_to_q calling dorgqr. LAPACK info = ", info
            endif
            Q = matrix(:k, :k)



        endfunction dgeq_to_q2


        function dgeq_to_q(A, tau) result(Q)
            real(dp), intent(in) :: A(:, :)
            real(dp), intent(in) :: tau(min(size(A, 1), size(A, 2)))

            integer  :: m, n, k
            real(dp) :: Q(size(A, 1), size(A, 1))

            integer  :: info
            real(dp) :: work(128*size(A, 1))
            integer  :: lwork

            m = size(A, 1); n = size(A, 2)
            k = size(tau)
            Q = id(m)
            info = 0
            work = 0.0_dp
            lwork = size(work)


            call dormqr('L', 'N', m, m, k, A, m, tau, Q, m, work, lwork, info)

        
        endfunction dgeq_to_q




        function jpvt_to_p(jpvt) result(P)
            integer :: jpvt(:)

            real(dp) :: P(size(jpvt), size(jpvt))
            integer  :: i

            P = 0.0_dp
            do i = 1, size(jpvt)
                P(jpvt(i), i) = 1.0_dp
            enddo


        endfunction jpvt_to_p






















        function gaussian() result(z)
            real(dp) :: z
            real(dp) :: theta, radius, r

            call random_number(r)
            theta = 2 * pi * r
            call random_number(r)
            radius = sqrt(-2.0_dp * log(r))
            ! z1 = radius * cos(theta)
            ! z2 = radius * sin(theta)
            z = radius * cos(theta)


        endfunction gaussian


        function gaussian_matrix(m, n) result(A)
            integer, intent(in) :: m
            integer, intent(in) :: n

            real(dp) :: A(m, n)
            integer  :: i, j

            do i = 1, m
                do j = 1, n
                    A(i, j) = gaussian()
                enddo
            enddo

            
        endfunction gaussian_matrix


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


        function id(n)
            integer, intent(in) :: n
            
            real(dp) :: id(n, n)
            integer  :: i

            id = 0.0_dp
            do i = 1, n
                id(i, i) = 1.0_dp
            enddo

        endfunction id


        subroutine print_matrix_nonzero(A, ounit)
            real(dp), intent(in) :: A(:, :)
            integer , intent(in) :: ounit
            integer :: i, j
            do i = 1, size(A, 1)
                do j = 1, size(A, 2)
                    if (abs(A(i, j)) .ge. 10e-4) then
                        write(ounit, "(i5, a, i5, f17.8)") i, " , ", j, A(i, j)
                    endif
                enddo
            enddo
        endsubroutine print_matrix_nonzero






endprogram main