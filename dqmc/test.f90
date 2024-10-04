program main
    use numbertypes
    implicit none

    integer , parameter :: n = 3
    real(dp)            :: A(n, n)
    real(dp)            :: B(n, n)
    integer             :: i, j

    A(:, 1) = [1, 2, 3]
    A(:, 2) = [4, 5, 6]
    A(:, 3) = [7, 8, 9]

    B(:, 1) = [-1, -2, -3]
    B(:, 2) = [-4, -5, -6]
    B(:, 3) = [-7, -8, -9]

    !do j = 1, n
    !    do i = 1, n
    !        A(i, j) = A(i, j) + B(j, i)
    !    enddo
    !enddo

    !do i = 1, n
    !    call daxpy(n, 1.0_dp, B(i, 1), n, A(1, i), 1)
    !enddo

    call daxpy(n*n, 1.0_dp, B, n, A, 1)


    print *, A(1, :)
    print *, A(2, :)
    print *, A(3, :)

    





endprogram main