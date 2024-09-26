module checkerboardtests_mod
    use numbertypes
    use customla_mod
    use checkerboard_mod
    implicit none

    contains


        subroutine right_matrixmultiply(A, B, n)
            !
            ! Updates:
            !
            ! A = A * B
            !
            ! where A and B are n x n matrices
            !
            real(dp), intent(inout) :: A
            real(dp), intent(in)    :: B
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
            real(dp), intent(inout) :: A
            real(dp), intent(in)    :: B
            integer , intent(in)    :: n

            real(dp), allocatable :: C(:, :)

            allocate(C(n, n))

            ! C = A
            call copy_matrix(A, C, n)
            ! A = B * C
            call dgemm('n', 'n', n, n, n, 1.0_dp, B, n, C, n, 0.0_dp, A, n)


        endsubroutine


        subroutine make_explicitcolor(color, A, n)
            !
            ! Sets A to the matrix represented by color explicitly.
            !
            ! color represents the matrix exponential of a strictly sparse matrix
            !
            type(ckbcolour), intent(in)  :: color
            real(dp)       , intent(out) :: A(n, n)
            integer        , intent(in)  :: n

            integer  :: k
            integer  :: i
            integer  :: j
            real(dp) :: d
            real(dp) :: ij
            real(dp) :: ji

            ! A = 0
            call zero_matrix(A, n)

            do k = 1, color%n
                i  = color%pairs(k)%i
                j  = color%pairs(k)%j
                ij = color%pairs(k)%ij
                ji = color%pairs(k)%ji
                d  = color%pairs(k)%d

                A(i, j) = ij
                A(j, i) = ji
                A(i, i) = d
                A(j, j) = d
            enddo


        endsubroutine make_explicitcolor


        subroutine make_explicitckb(ckb, A, n)
            !
            ! Sets A to the matrix exponential approximated by ckb (exactly the approximation)
            ! explicitly: by iterating through the sympairs stored in ckb.
            !
            ! ckb stores ckb%n colours ckb%colors(1), ..., ckb%colors(ckb%n)
            !
            ! A is made explicitly through the following process:
            !
            ! A = id
            ! A = A * ckb%colors(1)
            ! ...
            ! A = A * ckb%colors(ckb%n)
            !
            ! where the matrix exponential represented by ckb%colors(i)
            ! is made explicitly, then multiplied into A by dgemm
            !
            type(checkerboard), intent(in)    :: ckb
            real(dp)          , intent(out)   :: A(n, n)
            integer           , intent(in)    :: n

            integer :: i
            real(dp), allocatable :: B(:, :)

            allocate(B(n, n))

            ! A = id
            call make_identity(A, n)

            do i = 1, ckb%n
                ! B = ith color matrix exponential
                call make_explicitcolor(ckb%colours(i), B, n)
                ! A = A * B
                call right_matrixmultiply(A, B, n)
            enddo


        endsubroutine make_explicitckb


        subroutine make_explicitckbright(ckb, A, n)
            !
            ! Sets A to the matrix exponential approximated by ckb (exactly the approximation)
            ! by right multiplication by ckb: by calling the fast right ckb multiplication subroutine
            !
            type(checkerboard), intent(in)    :: ckb
            real(dp)          , intent(out)   :: A(n, n)
            integer           , intent(in)    :: n

            real(dp), allocatable :: work(:)

            allocate(work(n))

            ! A = id
            call make_identity(A, n)
            ! A = A * ckb
            call right_ckbmult(ckb, A, n, work)


        endsubroutine make_explicitckbright


        subroutine make_explicitckbleft(ckb, A, n)
            !
            ! Sets A to the matrix exponential approximated by ckb (exactly the approximation)
            ! by left multiplication by ckb: by calling the fast left ckb multiplication subroutine
            !
            type(checkerboard), intent(in)    :: ckb
            real(dp)          , intent(out)   :: A(n, n)
            integer           , intent(in)    :: n

            real(dp), allocatable :: work(:)

            allocate(work(n))

            ! A = id
            call make_identity(A, n)
            ! A = ckb * A
            call left_ckbmult(ckb, A, n, work)


        endsubroutine make_explicitckbleft



endmodule checkerboardtests_mod