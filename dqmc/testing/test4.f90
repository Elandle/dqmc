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
    integer         , parameter :: L        = 40
    real(dp)        , parameter :: dtau     = 0.08
    character(len=*), parameter :: ckbfname = "delta00.txt"



    type(checkerboard)    :: ckb
    integer               :: iounit, i, j, k
    real(dp)              :: T(n, n), A(n, n), B(n, n), G(n, n)
    real(dp)              :: VR(n, n), VL(n, n), wr(n), wi(n)
    real(dp)              :: beta


    beta = L * dtau
    ! T = hopping matrix
    call read_ckbT(T, n, ckbfname, iounit, 1.0_dp)
    ! A = exp(-beta * T)
    call expm(T, A, -beta)
    ! G = inv(id + exp(-beta * T))
    G = inverse(id(n) + A)
    call print_matrix(G, terminal)


    call diagonalize(T, VR, VL, wr, wi)
    
    write(terminal, "(a)") "T = "
    call print_matrix(T, terminal)
    write(terminal, "(a)") "Right eigenvectors = "
    call print_matrix(VR, terminal)
    write(terminal, "(a)") "Left eigenvectors = "
    call print_matrix(VL, terminal)
    write(terminal, "(a)") "Real part of eigenvalues = "
    call print_vector(wr, terminal)
    write(terminal, "(a)") "Imaginary part of eigenvalues = "
    call print_vector(wi, terminal)

    write(terminal, "(a)") "VR inner products:"
    do i = 1, n
        do j = 1, n
            write(terminal, "(i4, i4, f17.8)") i, j, dot_product(VR(:, i), VR(:, j))
        enddo
    enddo

    A = 0.0_dp
    do i = 1, n
        do j = 1, n
            do k = 1, n
                A(i, j) = A(i, j) + VR(i, k) * VR(j, k) * wr(k)
            enddo
        enddo
    enddo

    write(terminal, "(a)") "A = "
    call print_matrix(A, terminal)







    contains

        subroutine diagonalize(A, VR, VL, wr, wi)
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: VR(size(A, 1), size(A, 2))
            real(dp), intent(out) :: VL(size(A, 1), size(A, 2))
            real(dp), intent(out) :: wr(size(A, 1))
            real(dp), intent(out) :: wi(size(A, 1))

            real(dp), allocatable :: matrix(:, :)
            real(dp), allocatable :: work(:)
            integer               :: m, n, info, lwork

            m = size(A, 1); n = size(A, 2)
            lwork = 10 * m
            allocate(matrix(m, n))
            allocate(work(lwork))
            matrix = A


            call dgeev('V', 'V', m, matrix, m, wr, wi, VL, m, VR, m, work, lwork, info)

            deallocate(matrix)

        endsubroutine diagonalize



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