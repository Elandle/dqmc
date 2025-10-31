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

    

    integer, parameter :: n = 8
    type(checkerboard) :: ckb, ckbinv
    integer            :: iounit, i, j
    real(dp)           :: dtau
    real(dp)           :: A(n, n), expA(n, n), ckbA(n, n), R(n, n), work(n), ckbinvA(n, n), invexpA(n, n), Rinv(n, n)

    dtau = 1.0_dp
    call read_ckb(ckb, "delta02.txt", iounit, dtau)
    call read_ckb(ckbinv, "delta02.txt", iounit, -dtau)
    call read_ckbT(A, n, "delta02.txt", iounit, dtau)
    ckbA = 0.0_dp
    ckbinvA = 0.0_dp
    do i = 1, n
        ckbA(i, i) = 1.0_dp
        ckbinvA(i, i) = 1.0_dp
    enddo

    call right_ckbmult(ckb, ckbA, n, work)
    call right_ckbmult(ckbinv, ckbinvA, n, work)
    call expm(A, expA)
    invexpA = inverse(expA)
    write(terminal, "(a)") "A = "
    call print_matrix(A, terminal)
    write(terminal, "(a)") "ckbA = "
    call print_matrix(ckbA, terminal)
    write(terminal, "(a)") "expA = "
    call print_matrix(expA, terminal)
    write(terminal, "(a)") "ckbinvA = "
    call print_matrix(ckbinvA, terminal)
    write(terminal, "(a)") "invexpA = "
    call print_matrix(invexpA, terminal)
    R = expA - ckbA
    Rinv = invexpA - ckbinvA
    write(terminal, "(a)") "expA - ckbA = "
    call print_matrix(R, terminal)
    write(terminal, "(a)") "expA - ckbA norms:"
    write(terminal, "(a, f17.8)") "maxabs entry   = ", dmnorm(R, 'm')
    write(terminal, "(a, f17.8)") "one       norm = ", dmnorm(R, 'o')
    write(terminal, "(a, f17.8)") "infinity  norm = ", dmnorm(R, 'i')
    write(terminal, "(a, f17.8)") "frobenius norm = ", dmnorm(R, 'f')
    write(terminal, "(a, f17.8)") "two       norm = ", dmnorm(R, 't')
    write(terminal, "(a, f17.8)") "dtau   = ", dtau
    write(terminal, "(a, f17.8)") "dtau^2 = ", dtau * dtau
    write(terminal, "(a)") "invexpA - ckbinvA = "
    call print_matrix(Rinv, terminal)
    write(terminal, "(a)") "expA - ckbA norms:"
    write(terminal, "(a, f17.8)") "maxabs entry   = ", dmnorm(Rinv, 'm')
    write(terminal, "(a, f17.8)") "one       norm = ", dmnorm(Rinv, 'o')
    write(terminal, "(a, f17.8)") "infinity  norm = ", dmnorm(Rinv, 'i')
    write(terminal, "(a, f17.8)") "frobenius norm = ", dmnorm(Rinv, 'f')
    write(terminal, "(a, f17.8)") "two       norm = ", dmnorm(Rinv, 't')


    contains

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