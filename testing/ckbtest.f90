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

    

    integer         , parameter :: n        = 8
    real(dp)        , parameter :: dtau     = 0.08
    character(len=*), parameter :: ckbfname = "delta02.txt"



    type(checkerboard)    :: ckb
    integer               :: iounit, i, j, ii, jj
    real(dp)              :: d, ij, ji
    real(dp)              :: A(n, n), expA(n, n), ckbA(n, n), R(n, n), work(n)
    real(dp)              :: matrix(n, n)
    real(dp), allocatable :: colexps(:, :, :)


    call read_ckb(ckb, ckbfname, iounit, dtau)
    call read_ckbT(A, n, ckbfname, iounit, dtau)
    allocate(colexps(n, n, ckb%n))
    ckbA = 0.0_dp
    do i = 1, n
        ckbA(i, i) = 1.0_dp
    enddo

    write(terminal, "(a)") "A, its matrix exponential, and checkerboard exponential -------------"
    call right_ckbmult(ckb, ckbA, n, work)
    call expm(A, expA)
    write(terminal, "(a)") "A (dtau * read in matrix) = "
    call print_matrix(A, terminal)
    write(terminal, "(a)") "ckbA = "
    call print_matrix(ckbA, terminal)
    write(terminal, "(a)") "expA = "
    call print_matrix(expA, terminal)
    R = expA - ckbA
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
    write(terminal, "(a)") "----------------------------------------------------"


    write(terminal, "(a)") "Checkerboard ---------------------------------------"
    write(terminal, "(a, i4)") "ckb%n            = ", ckb%n
    do ii = 1, ckb%n
        write(terminal, "(a, i4, a, i4)") "ckb%colours(", ii, ")%n = ", ckb%colours(ii)%n
    enddo
    do ii = 1, ckb%n
        write(terminal, "(a, i4, a)") "ckb%colours(", ii, ") information ------------"
        do jj = 1, ckb%colours(ii)%n
            write(terminal, "(a, i4, a, i4, a)") "ckb%colours(", ii, ")%pairs(", jj, ") attributes:"
            i  = ckb%colours(1)%pairs(1)%i
            j  = ckb%colours(1)%pairs(1)%j
            d  = ckb%colours(1)%pairs(1)%d
            ij = ckb%colours(1)%pairs(1)%ij
            ji = ckb%colours(1)%pairs(1)%ji
            write(terminal, "(a, i4, a, i4)")                 "i, j      = ", i, ", ", j
            write(terminal, "(a, f17.8, a, f17.8, a, f17.8)") "d, ij, ji = ", d, ", ", ij, ", ", ji
            write(terminal, "(a, f17.8, a, f17.8)")           "            ", d , "  ", ij
            write(terminal, "(a, f17.8, a, f17.8)")           "exp       = ", ji, "  ", d
        enddo
        write(terminal, "(a, i4, a)") "ckb%colours(", ii, ") * id = "
        matrix = id(n)
        call left_colourmult(ckb%colours(ii), matrix, n, work)
        call print_matrix(matrix, terminal)
        write(terminal, "(a)") "nonzero entries:"
        call print_matrix_nonzero(matrix, terminal)
        write(terminal, "(a, i4, a)") "id * ckb%colours(", ii, ") = "
        matrix = id(n)
        call right_colourmult(ckb%colours(1), matrix, n, work)
        call print_matrix(matrix, terminal)
        write(terminal, "(a)") "nonzero entries:"
        call print_matrix_nonzero(matrix, terminal)
        colexps(:, :, ii) = matrix
    enddo
    write(terminal, "(a)") "colour1 * colour2 * ... * colourckb%n (dense multiplication) = "
    matrix = id(n)
    do i = 1, ckb%n
        matrix = matmul(matrix, colexps(:, :, i))
    enddo
    call print_matrix(matrix, terminal)
    write(terminal, "(a)") "ckbA - colour1 * ... = "
    call print_matrix(ckbA - matrix, terminal)





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