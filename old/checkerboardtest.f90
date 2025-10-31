program checkerboardtest
    use numbertypes
    use checkerboard_mod
    use checkerboardtests_mod
    use convenientla_mod
    implicit none

    type(checkerboard)    :: ckb, ckbinv
    character(len=100)    :: filename = "squareckb.txt"
    integer               :: iounit
    integer               :: n
    real(dp), allocatable :: A(:, :)
    real(dp), allocatable :: B(:, :)
    real(dp), allocatable :: C(:, :)
    real(dp), allocatable :: D(:, :)
    real(dp), allocatable :: E(:, :)
    integer i, j, k


    call read_ckb(ckb,    filename, iounit, 0.125_dp)
    call read_ckb(ckbinv, filename, iounit, 0.125_dp)
    n = 16

    allocate(A(n, n))
    allocate(B(n, n))
    allocate(C(n, n))
    allocate(D(n, n))
    allocate(E(n, n))

    call make_explicitckbright(ckb, A, n)
    call make_explicitckbleft(ckb, B, n)
    call make_explicitckb(ckb, C, n)

    print *, twonorm(A - B)
    print *, twonorm(C - A)
    print *, twonorm(B - C)

    call make_explicitckbright(ckbinv, D, n)

    call right_matrixmultiply(A, D, n)

    do j = 1, n
        do i = 1, n
            write(*, "(F12.6)", advance="no") A(i, j)
        enddo
        write(*, *) ""
    enddo

    








endprogram checkerboardtest