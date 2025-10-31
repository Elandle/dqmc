program main
    use iso_fortran_env, only: dp => real64, terminal => output_unit
    use checkerboard_mod
    use expm_mod
    implicit none

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

    

    integer, parameter   :: n = 8, LL = 5
    type(checkerboard)   :: ckb, ckbinv
    integer              :: iounit, i, j, l, lseed, h(n, LL), sigma
    real(dp)             :: dtau, alpha, mu, rand, U
    real(dp)             :: T(n, n), expT(n, n), ckbA(n, n), R(n, n), work(n), ckbinvA(n, n), invexpA(n, n), Rinv(n, n)
    real(dp)             :: Bckb(n, n), Bexact(n, n), matwork(n, n)
    integer, allocatable :: seed(:)
  
    call random_seed(size=lseed)
    allocate(seed(lseed))
    seed = 123
    call random_seed(put=seed)
    do i = 1, n
        do j = 1, LL
           call random_number(rand)
           if (rand .lt. 0.5_dp) then
              h(i, j) = -1
           else
              h(i, j) = 1
           endif
        enddo
    enddo
    write(terminal, "(a)") "h = "
    do i = 1, n
        write(terminal, "(*(i5))") h(i, :)
    enddo

    dtau  = 0.001_dp
    U     = 0.0_dp
    mu    = 0.0_dp
    alpha = acosh(exp(dtau * U / 2.0_dp))

    call read_ckb(ckb, "delta02.txt", iounit, dtau)
    call read_ckbT(T, n, "delta02.txt", iounit, dtau)
    call expm(T, expT)

    Bckb   = 0.0_dp
    Bexact = 0.0_dp
    do i = 1, n
        Bckb  (i, i) = 1.0_dp
        Bexact(i, i) = 1.0_dp
    enddo

    call random_number(Bckb)
    Bexact = Bckb


    l     = 1
    sigma = 1
    call right_Bmultckb(Bckb, ckb, work, l, h, sigma, alpha, dtau, mu)
    call right_Bmultexact(Bexact, expT, matwork, l, h, sigma, alpha, dtau, mu)

    write(terminal, "(a)") "Bckb = "
    call print_matrix(Bckb, terminal)
    write(terminal, "(a)") "Bexact = "
    call print_matrix(Bexact, terminal)


    contains

        subroutine right_Bmultckb(A, ckb, work, l, h, sigma, alpha, dtau, mu)
            real(dp)          , intent(inout) :: A(:, :)
            type(checkerboard), intent(in)    :: ckb
            real(dp)          , intent(out)   :: work(size(A, 1))
            integer           , intent(in)    :: l
            integer           , intent(in)    :: h(:, :)
            integer           , intent(in)    :: sigma
            real(dp)          , intent(in)    :: alpha
            real(dp)          , intent(in)    :: dtau
            real(dp)          , intent(in)    :: mu

            call right_diagmult(A, exp(sigma*alpha*h(:, l) + dtau*mu), size(A, 1))
            call right_ckbmult(ckb, A, size(A, 1), work)


        endsubroutine right_Bmultckb


        subroutine right_Bmultexact(A, expT, work, l, h, sigma, alpha, dtau, mu)
            real(dp), intent(inout) :: A(:, :)
            real(dp), intent(in)    :: expT(size(A, 1), size(A, 2))
            real(dp), intent(out)   :: work(size(A, 1), size(A, 2))
            integer , intent(in)    :: l
            integer , intent(in)    :: h(:, :)
            integer , intent(in)    :: sigma
            real(dp), intent(in)    :: alpha
            real(dp), intent(in)    :: dtau
            real(dp), intent(in)    :: mu

            call right_diagmult(A, exp(sigma*alpha*h(:, l)), size(A, 1))
            call right_matmul(A, expT, size(A, 1), work)
            ! TODO: combine with diagmult exp, once certain this is working
            A = exp(dtau * mu) * A


        endsubroutine right_Bmultexact
















        subroutine right_matmul(A, B, n, work)
            !
            ! Updates:
            !
            ! A = A * B
            !
            ! where A and B are n x n matrices.
            !
            ! Uses a supplied work matrix to hold a temporary copy of A (since
            ! there is no A = A * B general matrix update routine in BLAS/LAPACK).
            !
            ! This subroutine should be avoided at all costs, but it might be
            ! necessary to use at times.
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: B(n, n)
            integer , intent(in)    :: n
            real(dp), intent(in)    :: work(n, n)
            
            ! work = A
            call dlacpy('a', N, N, A, N, work, N)
            ! A = work * B
            call dgemm('n', 'n', N, N, N, 1.0_dp, work, N, B, N, 0.0_dp, A, N)


        endsubroutine right_matmul



        subroutine right_diagmult(A, D, n)
            !
            ! Updates:
            !
            ! A = A * D
            !
            ! Where D is a diagonal matrix stored as a vector
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: D(n)
            integer , intent(in)    :: n

            integer :: i

            ! TODO:
            ! Maybe change the do to a do concurrent?

            do i = 1, n
                call dscal(n, D(i), A(1, i), 1)
            enddo

            ! No BLAS:
            ! do i = 1, n
            !     A(:, i) = D(i) * A(:, i)
            ! enddo
            

        endsubroutine right_diagmult



        function inverse(A) result(invA)
            real(dp), intent(in)  :: A(:, :)

            real(dp) :: invA(size(A, 1), size(A, 2))
            integer  :: P(size(A, 1))
            real(dp) :: work(64 * size(A, 1))
            integer  :: info

            call dcopy(size(A, 1)*size(A, 1), A, 1, invA, 1)
            call dgetrf(size(invA, 1), size(invA, 1), invA, size(invA, 1), P, info)
            call dgetri(size(invA, 1), invA, size(invA, 1), P, work, 4*size(invA, 1), info)


        endfunction inverse

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




        subroutine print_matrix(A, ounit)
            real(dp), intent(in) :: A(:, :)
            integer , intent(in) :: ounit
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            do i = 1, n
                do j = 1, m
                    write(ounit, "(f17.8)", advance="no") A(i, j)
                enddo
                write(ounit, "(a)") ""
            enddo


        endsubroutine print_matrix




endprogram main