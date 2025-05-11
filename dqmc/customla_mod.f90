module customla_mod
    use numbertypes
    implicit none


    contains


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

        
        subroutine left_diagmult(A, D, n)
            !
            ! Updates:
            !
            ! A = D * A
            !
            ! Where D is a diagonal matrix stored as a vector
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: D(n)
            integer , intent(in)    :: n

            ! BEST
            ! For some reason, it seems only some BLAS/LAPACK distributions contain the
            ! dlascl2 subroutine

            ! call dlascl2(n, n, D, A, n)

            ! SECOND BEST
            ! Alternative BLAS:
            integer :: i
            do i = 1, n
                call dscal(n, D(i), A(i, 1), n)
            enddo

            ! THIRD BEST
            ! No BLAS:
            ! integer :: i
            ! do i = 1, n
            !     A(i, :) = D(i) * A(i, :)
            ! enddo
            
            
        endsubroutine left_diagmult


        subroutine right_diaginvmult(A, D, n)
            !
            ! Updates:
            !
            ! A = A * inv(D)
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
                call dscal(n, 1.0_dp / D(i), A(1, i), 1)
            enddo

            ! No BLAS:
            ! do i = 1, n
            !     A(:, i) = A(:, i) / D(i)
            ! enddo
            

        endsubroutine right_diaginvmult


        subroutine left_diaginvmult(A, D, n)
            !
            ! Updates:
            !
            ! A = inv(D) * A
            !
            ! Where D is a diagonal matrix stored as a vector
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: D(n)
            integer , intent(in)    :: n

            ! For some reason, it seems only some BLAS/LAPACK distributions contain the
            ! dlarscl2 subroutine

            ! BEST
            ! call dlarscl2(n, n, D, A, n)

            ! SECOND BEST
            ! Alternative BLAS:
            integer :: i
            do i = 1, n
                call dscal(n, 1.0_dp / D(i), A(i, 1), n)
            enddo

            ! THIRD BEST
            ! No BLAS:
            ! integer :: i
            ! do i = 1, n
            !     A(i, :) = A(i, :) / D(i)
            ! enddo
            
            
        endsubroutine left_diaginvmult


        subroutine diag(A, D, n)
        !
        ! Sets:
        !
        ! D = diag(A)
        !
        ! Where A is an n x n matrix and D is stored as a length n vector
        !
        real(dp), intent(in)  :: A(n, n)
        real(dp), intent(out) :: D(n)
        integer , intent(in)  :: n

        integer :: i

        ! TODO:
        ! Test the BLAS version

        do concurrent (i = 1 : n)
            D(i) = A(i, i)
        enddo

        ! BLAS version:
        ! call dcopy(n, A(1, 1), n+1, D, 1)


        endsubroutine diag


        subroutine uppertri(A, B, n)
            !
            ! Sets:
            !
            ! B = uppertri(A)
            !
            ! That is, B is set to be the upper triangular part of A (including the diagonal),
            ! with all other entries set to zero.
            !
            real(dp), intent(in)  :: A(n, n)
            real(dp), intent(out) :: B(n, n)
            integer , intent(in)  :: n

            ! Zero the lower triangular part of B (including the diagonal)
            call dlaset('L', N, N, 0.0_dp, 0.0_dp, B, N)

            ! Copy the upper triangular part of A (including the diagonal) to B
            call dlacpy('U', N, N, A, N, B, N)


        endsubroutine uppertri


        ! TODO:
        ! Compare the two permutation routines
        ! Get rid of gotos in second permutation routine
        ! Better routine?


        subroutine invert_permutation_old(P, I, n)
            !
            ! Sets:
            !
            ! I = inv(P)
            !
            ! In terms of permutation matrices stored as length n vectors
            !
            integer, intent(in)  :: P(n)
            integer, intent(out) :: I(n)
            integer, intent(in)  :: n

            integer :: j

            do j = 1, n
                I(P(j)) = j ! Better way to do this using array indexing?
            enddo

            ! Maybe:
            ! I(P) = [(j, j = 1, n)]


        endsubroutine invert_permutation_old


        subroutine invert_permutation(P, n)
            ! Knuth, The Art of Computer Programming, Volume 1
            ! Section 1.3.3, algorithm I
            integer, intent(inout) :: P(n)
            integer, intent(in)    :: n

            integer :: m, j, i

            ! I1
            m = n
            j = -1

            ! I2
20          i = P(m)
            if (i .lt. 0) then
                goto 50
            endif

            ! I3
30          P(m) = j
            j = -m
            m = i
            i = P(m)

            ! I4
            if (i .gt. 0) then
                goto 30
            else
                i = j
            endif

            ! I5
50          P(m) = -i

            ! I6
            m = m - 1
            if (m .gt. 0) then
                goto 20
            endif


        endsubroutine invert_permutation













        subroutine permutecols(A, P, n)
            !
            ! Updates:
            !
            ! A = A * P
            !
            ! Where P is a permutation matrix stored as a vector.
            !
            ! This subroutine is implemented basically as a copy of Julia's intrinsic permutecols
            !
            real(dp), intent(inout) :: A(n, n)
            integer , intent(inout) :: P(n)
            integer , intent(in)    :: n

            integer :: count
            integer :: start
            integer :: ptr
            integer :: next

            count = 0; start = 0
        
            do while (count .lt. n)
                start = findloc(P, start + 1, mask = P .ne. 0, dim = 1)
                if (start .eq. 0) exit
                ptr = start

                ! Equivalent (without using findloc):
                ! do start = start + 1, n
                !     if (P(start) .ne. 0) then
                !         ptr = start
                !         exit
                !     endif
                ! enddo
        
                next = P(start)
                count = count + 1
        
                do while (next .ne. start)
                    call dswap(n, A(1, ptr), 1, A(1, next), 1) ! Swap columns A(:, ptr) <--> A(:, next)
                    P(ptr) = 0
                    ptr = next
                    next = P(next)
                    count = count + 1
                end do

                P(ptr) = 0
            end do


        endsubroutine permutecols


        ! These following permutation algorithms are all based off the same underlying algorithm
        ! TODO:
        ! Find an inplace permutation algorithm that doesn't need an inverse permutation


        subroutine permute_matrix_columns(A, m, P, Pinv)
            real(dp), intent(inout) :: A(m, m)
            integer , intent(in)    :: m
            integer , intent(inout) :: P(m)
            integer , intent(inout) :: Pinv(m)

            integer :: i, j, k

            do i = 1, m
                j = P(i)
                if (j .ne. i) then
                    k = Pinv(i)
                    do
                        if (i .lt. j .and. i .lt. k) then
                            if (j .eq. k) then
                                call rotate_matrix_columns(P, i, A, m)
                                exit
                            endif
                            j = P(j)
                            if (j .eq. k) then
                                call rotate_matrix_columns(P, i, A, m)
                                exit
                            endif
                            k = Pinv(k)
                        else
                            exit
                        endif
                    enddo
                endif
            enddo


        endsubroutine permute_matrix_columns

        subroutine permute_matrix_rows(A, m, P, Pinv)
            real(dp), intent(inout) :: A(m, m)
            integer , intent(in)    :: m
            integer , intent(inout) :: P(m)
            integer , intent(inout) :: Pinv(m)

            integer :: i, j, k

            do i = 1, m
                j = P(i)
                if (j .ne. i) then
                    k = Pinv(i)
                    do
                        if (i .lt. j .and. i .lt. k) then
                            if (j .eq. k) then
                                call rotate_matrix_rows(P, i, A, m)
                                exit
                            endif
                            j = P(j)
                            if (j .eq. k) then
                                call rotate_matrix_rows(P, i, A, m)
                                exit
                            endif
                            k = Pinv(k)
                        else
                            exit
                        endif
                    enddo
                endif
            enddo


        endsubroutine permute_matrix_rows

        subroutine rotate_matrix_columns(P, leader, A, m)
            ! Part of permute
            integer , intent(inout) :: P(m)
            integer , intent(in)    :: leader
            real(dp), intent(inout) :: A(m, m)
            integer , intent(in)    :: m

            integer :: i

            i = P(leader)
            do
                if (i .eq. leader) then
                    exit
                else
                    call dswap(m, A(1, i), 1, A(1, leader), 1)
                    i = P(i)
                endif
            enddo

        endsubroutine rotate_matrix_columns


        subroutine rotate_matrix_rows(P, leader, A, m)
            ! Part of permute
            integer , intent(inout) :: P(m)
            integer , intent(in)    :: leader
            real(dp), intent(inout) :: A(m, m)
            integer , intent(in)    :: m

            integer :: i

            i = P(leader)
            do
                if (i .eq. leader) then
                    exit
                else
                    call dswap(m, A(i, 1), m, A(leader, 1), m)
                    i = P(i)
                endif
            enddo

        endsubroutine rotate_matrix_rows




        subroutine permute(P, Pinv, x, m)
            ! Fich, Munro, and Poblete
            ! Permuting in Place
            ! Figure 5
            integer , intent(inout) :: P(m)
            integer , intent(inout) :: Pinv(m)
            real(dp), intent(inout) :: x(m)
            integer , intent(in)    :: m

            integer :: i, j, k

            do i = 1, m
                j = P(i)
                if (j .ne. i) then
                    k = Pinv(i)
                    do
                        if (i .lt. j .and. i .lt. k) then
                            if (j .eq. k) then
                                call rotate(P, i, x)
                                exit
                            endif
                            j = P(j)
                            if (j .eq. k) then
                                call rotate(P, i, x)
                                exit
                            endif
                            k = Pinv(k)
                        else
                            exit
                        endif
                    enddo
                endif
            enddo


        endsubroutine permute

        subroutine rotate(P, leader, x)
            ! Part of permute
            integer , intent(inout) :: P(:)
            integer , intent(in)    :: leader
            real(dp), intent(inout) :: x(:)

            integer :: i

            i = P(leader)
            do
                if (i .eq. leader) then
                    exit
                else
                    call swap(x, i, leader)
                    i = P(i)
                endif
            enddo

        endsubroutine rotate

        subroutine swap(x, i, j)
            real(dp), intent(inout) :: x(:)
            integer , intent(in)    :: i
            integer , intent(in)    :: j

            real(dp) :: temp

            temp = x(i)
            x(i) = x(j)
            x(j) = temp

        endsubroutine swap







        subroutine invert(A, n, P, work, lwork, info)
            !
            ! Sets:
            !
            ! A = inv(A)
            !
            ! where A is a general matrix.
            !
            ! A is inverted in two steps:
            !
            ! First, A is A = P * L * U factorised by dgetrf
            ! Then, A is inverted by dgetri (using the result of dgetrf)
            !
            real(dp), intent(inout) :: A(n, n)
            integer , intent(in)    :: n
            integer , intent(out)   :: P(n)
            real(dp), intent(out)   :: work(lwork)
            integer , intent(in)    :: lwork
            integer , intent(out)   :: info

            ! A = P * L * U
            call dgetrf(n, n, A, n, P, info)
            ! A = inv(A)
            call dgetri(n, A, n, P, work, lwork, info)


        endsubroutine invert

        
        subroutine make_identity(A, n)
            !
            ! Sets:
            !
            ! A = id
            !
            ! where id is the identity matrix and A is an n x n matrix
            !
            real(dp), intent(inout) :: A(n, n)
            integer , intent(in)    :: n

            ! A = id
            call dlaset('A', n, n, 0.0_dp, 1.0_dp, A, n)
            

        endsubroutine make_identity


        subroutine zero_matrix(A, n)
            !
            ! Sets:
            !
            ! A = 0
            !
            real(dp), intent(out) :: A(n, n)
            integer , intent(in)  :: n

            ! A = 0
            call dlaset('A', n, n, 0.0_dp, 0.0_dp, A, n)


        endsubroutine zero_matrix


        subroutine copy_matrix(A, B, n)
            !
            ! Sets:
            !
            ! B = A
            !
            real(dp), intent(in)  :: A(n, n)
            real(dp), intent(out) :: B(n, n)
            integer , intent(in)  :: n

            ! B = A
            call dlacpy('A', n, n, A, n, B, n)


        endsubroutine copy_matrix


        subroutine add_trans(A, B, n)
            !
            ! Updates:
            !
            ! A = A + trans(B)
            !
            ! where A and B are n x n matrices.
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: B(n, n)
            integer , intent(in)    :: n

            integer :: i

            ! Intel's MKL might have its own version of this operation

            ! BLAS:
            do i = 1, n
                call daxpy(n, 1.0_dp, B(i, 1), n, A(1, i), 1)
            enddo

            ! No BLAS:
            ! do i = 1, n
            !     A(:, i) = A(:, i) + B(i, :)
            ! enddo


        endsubroutine add_trans


        subroutine add_matrix(A, B, n)
            !
            ! Updates:
            !
            ! A = A + B
            !
            ! where A and B are n x n matrices
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: B(n, n)
            integer , intent(in)    :: n

            ! BLAS:
            call daxpy(n*n, 1.0_dp, B, 1, A, 1)

            ! No BLAS:
            ! integer :: i, j
            ! do j = 1, n
            !     do i = 1, n
            !         A(i, j) = A(i, j) + B(i, j)
            !     enddo
            ! enddo


        endsubroutine add_matrix


        subroutine left_matmul(A, B, n, work)
            !
            ! Updates:
            !
            ! A = B * A
            !
            ! where A and B are n x n matrices.
            !
            ! Uses a supplied work matrix to hold a temporary copy of A (since
            ! there is no A = B * A general matrix update routine in BLAS/LAPACK).
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
            ! A = B * work
            call dgemm('n', 'n', N, N, N, 1.0_dp, B, N, work, N, 0.0_dp, A, N)


        endsubroutine left_matmul


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


        subroutine trans(A, B, n)
            !
            ! Sets:
            !
            ! A = trans(B)
            !
            real(dp), intent(out) :: A(n, n)
            real(dp), intent(in)  :: B(n, n)
            integer , intent(in)  :: n

            integer :: i
            integer :: j

            ! TODO: see if BLAS/LAPACK can be implemented

            do j = 1, n
                do i = 1, n
                    A(i, j) = B(j, i)
                enddo
            enddo


        endsubroutine trans


        subroutine dlaswpc(n, a, lda, k1, k2, ipiv, incx)
            !
            !    MODIFIED VERSION OF DLASWP
            !  -- LAPACK auxiliary routine --
            !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
            !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
            !
            !     .. Scalar Arguments ..
            integer            incx, k1, k2, lda, n
            !     ..
            !     .. Array Arguments ..
            integer            ipiv(*)
            double precision   A(lda, *)
            !     ..
            !
            ! =====================================================================
            !
            !     .. Local Scalars ..
            integer            i, i1, i2, inc, ip, ix, ix0, j, k, n32
            double precision   temp
            !     ..
            !     .. Executable Statements ..
            !
            !     Interchange row i with row ipiv(k1+(i-k1)*abs(incx)) for each of rows
            !     k1 through k2.
            !
            if (incx .gt. 0) then
                ix0 = k1
                i1  = k1
                i2  = k2
                inc = 1
            elseif (incx .lt. 0) then
                ix0 = k1 + (k1-k2) * incx
                i1  = k2
                i2  = k1
                inc = -1
            else
                return
            endif
            !
            n32 = (n/32) * 32
            if (n32 .ne. 0) then
                do 30 j = 1, n32, 32
                    ix = ix0
                    do 20 i = i1, i2, inc
                        ip = ipiv(ix)
                        if (ip .ne. i) then
                            do 10 k = j, j + 31
                                temp     = a(k , i)
                                a(k , i) = a(k, ip)
                                a(k, ip) = temp
        10                  continue
                        endif
                        ix = ix + incx
        20          continue
        30      continue
            endif
        
            if (n32 .ne. n) then
                n32 = n32 + 1
                ix  = ix0
                do 50 i = i1, i2, inc
                    ip = ipiv(ix)
                    if (ip .ne. i) then
                        do 40 k = n32, n
                            temp     = a(k , i)
                            a(k , i) = a(k, ip)
                            a(k, ip) = temp
        40              continue
                    endif
                    ix = ix + incx
        50      continue
            endif
            !
            return
            !
            !     End of DLASWP
            !


        endsubroutine dlaswpc


        subroutine colpivswap(A, piv, n, matwork)
            !
            ! Swaps the columns of the n x n matrix A according to the n
            ! long integer vector piv
            !
            ! If piv(i) = k, then column i of A becomes column k of A
            !
            ! In other words, for i = 1, 2, ..., n:
            !
            !       A(:, piv(i)) = A(:, i)
            !
            ! Where this replacement is done independently (changes where
            ! columns are do not affect other column changes)
            !
            real(dp), intent(inout) :: A(n, n)
            integer , intent(in)    :: piv(n)
            integer , intent(in)    :: n
            real(dp), intent(out)   :: matwork(n, n)

            integer i

            !
            ! TODO:
            ! Can this be done without copying A?
            !

            call copy_matrix(A, matwork, n)

            do i = 1, n
               if (piv(i) .ne. i) then ! Do not copy a column if it is already in the right place
                   ! A(:, piv(i)) = matwork(:, i)
                   call dcopy(n, matwork(1, i), 1, A(1, piv(i)), 1)
               endif
            enddo


        endsubroutine colpivswap


        function matdiff(A, B, n)
            real(dp), intent(in) :: A(n, n)
            real(dp), intent(in) :: B(n, n)
            integer , intent(in) :: n
            real(dp)             :: matdiff
            real(dp), external   :: dnrm2

            matdiff = dnrm2(n*n, A - B, 1) / n


        endfunction matdiff


        subroutine determinant(det, A, n, P, info)
            !
            ! Sets:
            !
            ! det = det(A)
            !
            ! by first factorising A = P * L * U (destroying A),
            ! and calculating the determinant by multiplying the diagonal of U
            ! and adjusting the sign based on P
            !
            real(dp)                :: det
            real(dp), intent(inout) :: A(n, n)
            integer , intent(in)    :: n
            integer , intent(out)   :: P(n)
            integer , intent(out)   :: info

            integer :: i

            ! A = P * L * U
            call dgetrf(n, n, A, n, P, info)

            det = A(1, 1)
            if (P(1) .ne. 1) then
                det = -1 * det
            endif

            do i = 2, n
                det = det * A(i, i)
                if (P(i) .ne. i) then
                    det = -1 * det
                endif
            enddo


        endsubroutine determinant


        subroutine zdeterminant(det, A, n, P, info)
            !
            ! Sets:
            !
            ! det = det(A)
            !
            ! by first factorising A = P * L * U (destroying A),
            ! and calculating the determinant by multiplying the diagonal of U
            ! and adjusting the sign based on P
            !
            complex(dp)                :: det
            complex(dp), intent(inout) :: A(n, n)
            integer    , intent(in)    :: n
            integer    , intent(out)   :: P(n)
            integer    , intent(out)   :: info

            integer :: i

            ! A = P * L * U
            call zgetrf(n, n, A, n, P, info)

            det = A(1, 1)
            if (P(1) .ne. 1) then
                det = -1 * det
            endif

            do i = 2, n
                det = det * A(i, i)
                if (P(i) .ne. i) then
                    det = -1 * det
                endif
            enddo


        endsubroutine zdeterminant


    subroutine eigenvalues(A, wr, wi)
        real(dp), intent(in)  :: A(:, :)
        real(dp), intent(out) :: wr(size(A, 1))
        real(dp), intent(out) :: wi(size(A, 1))

        integer :: m
        real(dp), allocatable :: B(:, :)
        real(dp), allocatable :: work(:)
        integer               :: lwork
        integer               :: info

        m = size(A, 1)
        lwork = 8 * m
        allocate(B(m, m))
        allocate(work(lwork))
        call copy_matrix(A, B, m)


        call dgeev('N', 'N', m, B, m, wr, wi, A, m, A, m, work, lwork, info)


    endsubroutine eigenvalues

endmodule customla_mod