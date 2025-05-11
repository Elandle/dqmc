program main
    use iso_fortran_env, only: dp => real64, terminal => output_unit
    implicit none
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

    interface
        function dlange(norm, m, n, A, lda, work) result(res)
            import :: dp
            character(len=*), intent(in) :: norm
            integer, intent(in) :: m, n, lda
            real(dp), intent(in) :: A(lda, *)
            real(dp), intent(out) :: work(*)
            real(dp) :: res
        end function dlange
    end interface


    integer :: m, n, k, pp
    real(dp), allocatable :: A(:, :), P(:, :), Q(:, :), R(:, :)

    k  = 600
    pp = 8

    call readin_mtx("ex33.mtx", A, m, n)
    allocate(P(n, n))
    allocate(Q(n, n))
    allocate(R(k, n))


    call ssrqrcp(A, k, pp, P, Q, R)

    ! call print_matrix(matmul(A, P), terminal, "AP = ")
    ! call print_matrix(matmul(Q(:, :k), R), terminal, "Q(:, :k)R = ")
    ! write(terminal, "(a, f17.8)") "frobenius-norm relative error of approximation = ", &
    ! & dmnorm(matmul(A, P) - matmul(Q(:, :k), R), "f") / dmnorm(matmul(A, P), "f")
    write(terminal, "(a)") "k, relative error"
    write(terminal, "(i4, f17.8)") k, dmnorm(matmul(A, P) - matmul(Q(:, :k), R), "f") / dmnorm(matmul(A, P), "f")



    contains


        subroutine ssrqrcp(A, k, pp, P, Q, R)
            !
            ! A: m x n
            ! P: n x n
            ! Q: m x m
            ! R: k x n
            !
            ! l = k + pp
            ! omega      : l x m
            ! B          : l x n
            ! A1         : m x n
            ! Qb         : l x l
            ! Rb         : l x n
            ! A1(:, :k)  : m x k
            ! Q(:, :k)   : m x k --> trans(Q(:, :k)): k x m
            ! A1(:, k+1:): m x (n-k)
            ! R11        : m x k
            ! R12        : k x (n-k)
            !
            ! R(:, :k)   = R11(:k, :)
            ! R(:, k+1:) = R12
            !
            ! subroutine qr(A, Q, R)
                !
                ! A: m x n
                ! Q: m x m
                ! R: m x n
                !
            ! subroutine qrcp(A, P, Q, R)
                !
                ! A: m x n
                ! P: n x n
                ! Q: m x m
                ! R: m x n
                !
            real(dp), intent(in)    :: A(:, :)
            integer , intent(in)    :: k
            integer , intent(in)    :: pp
            real(dp), intent(inout) :: P(size(A, 2), size(A, 2))
            real(dp), intent(out)   :: Q(size(A, 1), size(A, 1))
            real(dp), intent(out)   :: R(k, size(A, 2))

            integer  :: m, n, l
            real(dp) :: omega(k+pp, size(A, 1))
            real(dp) :: B(k+pp, size(A, 2))
            real(dp) :: Qb(k+pp, k+pp), Rb(k+pp, size(A, 2))
            real(dp) :: R11(size(A, 1), k), R12(k, size(A, 2)-k)
            real(dp) :: A0(size(A, 1), size(A, 2))
            real(dp) :: A1(size(A, 1), size(A, 2))

            m = size(A, 1); n = size(A, 2)
            l = size(omega, 1)
            A0 = A

            
            omega = gaussian_matrix(l, m)
            B = matmul(omega, A0)
            call qrcp(B, P, Qb, Rb)
            A1 = matmul(A0, P)
            call qr(A1(:, :k), Q, R11)
            R12 = matmul(transpose(Q(:, :k)), A1(:, k+1:))
            R(:, :k)   = R11(:k, :)
            R(:, k+1:) = R12


        endsubroutine ssrqrcp


        subroutine qr(A, Q, R)
            !
            ! A: m x n
            ! Q: m x m
            ! R: m x n
            !
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: Q(size(A, 1), size(A, 1))
            real(dp), intent(out) :: R(size(A, 1), size(A, 2))

            real(dp) :: matrix(size(A, 1), size(A, 2))
            integer  :: m, n

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
            Q(:, :) = dgeq_to_q(matrix, tau)
            R(:, :) = dgeq_to_r(matrix)

            
        endsubroutine qr


        subroutine qrcp(A, P, Q, R)
            !
            ! A: m x n
            ! P: n x n
            ! Q: m x m
            ! R: m x n
            !
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: P(size(A, 2), size(A, 2))
            real(dp), intent(out) :: Q(size(A, 1), size(A, 1))
            real(dp), intent(out) :: R(size(A, 1), size(A, 2))

            real(dp) :: matrix(size(A, 1), size(A, 2))
            integer  :: m, n

            integer  :: jpvt(size(A, 2))
            real(dp) :: tau(min(size(A, 1), size(A, 2)))
            real(dp) :: work(10*size(A, 2))
            integer  :: lwork
            integer  :: info

            m = size(A, 1); n = size(A, 2)
            lwork = size(work)
            info = 0
            jpvt = 0
            tau = 0.0_dp
            work = 0.0_dp

            matrix = A

            call dgeqp3(m, n, matrix, m, jpvt, tau, work, lwork, info)
            if (info .lt. 0) then
                write(terminal, "(a, i4)") "Error in subroutine qrcp calling dgeqp3. LAPACK info = ", info
            endif
            Q = dgeq_to_q(matrix, tau)
            R = dgeq_to_r(matrix)
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
            Q(:, :) = id(m)
            info = 0
            work = 0.0_dp
            lwork = size(work)


            call dormqr('L', 'N', m, m, k, A, m, tau, Q, m, work, lwork, info)

        
        endfunction dgeq_to_q




        function jpvt_to_p(jpvt) result(P)
            integer :: jpvt(:)

            real(dp) :: P(size(jpvt), size(jpvt))
            integer  :: i

            P(:, :) = 0.0_dp
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


        real(dp) function dmnorm(A, norm)
        !
        ! max(abs(A(i,j))), NORM = 'M' or 'm'
        ! norm1(A),         NORM = '1', 'O' or 'o'
        ! normI(A),         NORM = 'I' or 'i'
        ! normF(A),         NORM = 'F', 'f', 'E' or 'e'
        ! norm2(A),         NORM = 'T', 't', or '2'
        ! where  norm1  denotes the  one norm of a matrix (maximum column sum),
        ! normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
        ! normF  denotes the  Frobenius norm of a matrix (square root of sum of
        ! squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
        !
        real(dp) , intent(in) :: A(:, :)
        character, intent(in) :: norm

        real(dp), allocatable :: work(:)
        integer :: m, n, i, j

        m = size(A, 1)
        n = size(A, 2)

        if ((norm .eq. 'I') .or. (norm .eq. 'i')) then
            allocate(work(m))
        endif

        if ((norm .eq. '2') .or. (norm .eq. 'T') .or. (norm .eq. 't')) then
            dmnorm = twonorm(A)
        else
            dmnorm = dlange(norm, m, n, A, m, work)
        endif

        if (allocated(work)) then
            deallocate(work)
        endif

    endfunction dmnorm

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

        deallocate(B)
        deallocate(S)
        deallocate(work)

        
    endfunction twonorm


    subroutine readin_mtx(filename, A, m, n)
        implicit none
        character(len=*), intent(in) :: filename
        real(8), allocatable, intent(out) :: A(:,:)
        integer, intent(out) :: m, n
    
        integer :: L, iunit, i, row_idx, col_idx, ios
        real(8) :: val
        character(len=256) :: line
    
        iunit = 10
        open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
    

        do
            read(iunit, '(A)', iostat=ios) line
            if (ios .ne. 0) stop
            if (line(1:1) .ne. '%') exit
        end do
    

        read(line, *) m, n, L
    
        allocate(A(m, n))
        A = 0.0_dp
    
        ! Read and fill matrix
        do i = 1, L
            read(iunit, *) row_idx, col_idx, val
            A(row_idx, col_idx) = val
        end do
    
        close(iunit)
    end subroutine readin_mtx
    
    





endprogram main