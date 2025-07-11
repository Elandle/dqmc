program main
    use numbertypes
    use simulate_mod
    use printing_mod
    use iso_fortran_env, only: stdout => output_unit, stdin => input_unit
    implicit none
    type(Simulation)      :: S
    integer               :: iunit
    integer               :: i, l, j, k
    real(dp), allocatable :: A(:, :), B(:, :), C(:, :), D(:, :), id(:, :)
    real(dp)              :: R

    call setup_simulation_input(S, "input.txt", iunit, output_unit)
    S%h = 1
    ! Open debugging output file
    open(newunit=S%dunit, file=S%debfilename, action="write", status="replace")

    write(stdout, "(a)") "Calling simulate..."
    call simulate(S)
    write(stdout, "(a)") "simulate done."

    allocate(A(S%N, S%N)); allocate(B(S%N, S%N)); allocate(C(S%N, S%N))
    allocate(D(S%N, S%N)); allocate(id(S%N, S%N))

    id = 0.0_dp
    do j = 1, S%N
        id(j, j) = 1.0_dp
    enddo

    call print_integer_matrix(S%h, stdout, "Hubbard Stratonovich field h = ")
    i = 2; l = 3
    write(stdout, "(a, i4, a, i4, a)") "Flipping h(i, l) = h(", i, ", ", l, ")..."
    S%h(i, l) = -S%h(i, l)
    call print_integer_matrix(S%h, stdout, "Hubbard Stratonovich field h = ")
    write(stdout, "(a)") "Constructing new Gup..."
    call newG(S, l, 1)
    call print_matrix(S%Gup, stdout, "Before flip update, Gup = ")
    ! call flipupdate5(S, i, l, 1)
    call print_matrix(S%Gup, stdout, "After flip update, Gup = ")
    A = S%Gup
    write(stdout, "(a)") "Constructing new Gup..."
    call newG(S, l, 1)
    call print_matrix(S%Gup, stdout, "New Gup = ")
    call print_matrix(A - S%Gup, stdout, "Flip updated - new Gup = ")





    ! Agreement here
    ! S%h(i, l) = -S%h(i, l) ! undo flip
    ! call newG(S, l, 1)     ! Unflipped G
    ! D = 0.0_dp
    ! D(i, i) = exp(-2 * 1 * S%alpha * S%h(i, l)) - 1.0_dp
    ! A = matmul(S%Gup, inverse(id + matmul(D, id-S%Gup)))
    ! S%h(i, l) = -S%h(i, l)
    ! call newG(S, l, 1)
    ! call print_matrix(A - S%Gup, stdout, "updated - new Gup = ")

    ! Agreement here
    ! S%h(i, l) = -S%h(i, l) ! undo flip
    ! call newG(S, l, 1)     ! Unflipped G
    ! D = 0.0_dp
    ! D(i, i) = exp(-2 * 1 * S%alpha * S%h(i, l)) - 1.0_dp
    ! A = id + matmul(D, id-S%Gup)
    ! R = 1.0_dp + D(i, i) * (1.0_dp - S%Gup(i, i))
    ! B = id - (1.0_dp/R) * matmul(D, (id - S%Gup))
    ! AB should be id
    ! call print_matrix(matmul(A, B), stdout, "AB = ")


    ! Agreement
    ! S%h(i, l) = -S%h(i, l) ! undo flip
    ! call newG(S, l, 1)     ! Unflipped G
    ! D = 0.0_dp
    ! D(i, i) = exp(-2 * 1 * S%alpha * S%h(i, l)) - 1.0_dp
    ! R = 1.0_dp + D(i, i) * (1.0_dp - S%Gup(i, i))
    ! A = S%Gup - (1.0_dp/R) * matmul(S%Gup, matmul(D, id - S%Gup))
    ! S%h(i, l) = -S%h(i, l)
    ! call newG(S, l, 1)
    ! call print_matrix(A - S%Gup, stdout, "updated - new Gup = ")



    S%h(i, l) = -S%h(i, l) ! undo flip

    i = 1; l = 1

    call newG(S, l, 1)     ! Unflipped G
    D = 0.0_dp
    D(i, i) = exp(-2 * 1 * S%alpha * S%h(i, l)) - 1.0_dp
    R = 1.0_dp + D(i, i) * (1.0_dp - S%Gup(i, i))
    C = id - S%Gup
    do j = 1, S%N
        do k = 1, S%N
            A(j, k) = S%Gup(j, k) - (1.0_dp/R) * S%Gup(j, i) * D(i, i) * C(i, k)
        enddo
    enddo
    S%h(i, l) = -S%h(i, l)
    call newG(S, l, 1)
    call print_matrix(A, stdout, "flip updated Gup = ")
    call print_matrix(S%Gup, stdout, "New Gup = ")
    call print_matrix(A - S%Gup, stdout, "updated - new Gup = ")




    call print_matrix(S%T, stdout, "Hopping matrix T = ")
    ! Agreement on expT
    call print_matrix(S%expT, stdout, "exp(dtau*T) = ")
    write(stdout, "(a, f17.8)") "dtau = ", S%dtau
    write(stdout, "(a, f17.8)") "U = ", S%U
    write(stdout, "(a, f17.8)") "mu = ", S%mu
    ! Agreement on alpha
    write(stdout, "(a, f17.8)") "alpha = arccosh(exp(dtau * U / 2)) = ", S%alpha

    ! Agreement here
    ! do l = 1, S%L
        ! call make_B(S, C, l, 1)
        ! write(stdout, "(a, i4)") "make_B (up) l = ", l
        ! call print_matrix(C, stdout)

        ! call make_A(S, A, l, 1)
        ! B = exp(S%dtau * S%mu) * matmul(A, S%expT)
        ! write(stdout, "(a, i4)") "main make_B (up) l = ", l
        ! call print_matrix(B, stdout)

        ! call print_matrix(B - C, stdout, "main - make_B = ")
    ! enddo

    ! Agreement here
    ! newG test for:
    ! L = 4
    ! l = 2
    ! G(2) = inv(id + B(2)B(1)B(4)B(3))
    ! call make_B(S, A, 2, 1) ! A = B(2)
    ! call make_B(S, B, 1, 1) ! B = B(1)
    ! A = matmul(A, B)        ! A = AB = B(2)B(1)
    ! call make_B(S, B, 4, 1) ! B = B(4)
    ! A = matmul(A, B)        ! A = AB = B(2)B(1)B(4)
    ! call make_B(S, B, 3, 1) ! B = B(3)
    ! A = matmul(A, B)        ! A = AB = B(2)B(1)B(4)B(3)
    ! A = id + A
    ! do i = 1, S%N
        ! A(i, i) = A(i, i) + 1.0_dp
    ! enddo
    ! A = inv(A) = inv(id + B(2)B(1)B(4)B(3)) = G(2)
    ! A = inverse(A)
    ! call print_matrix(A, stdout, "main computed G(2) = ")
    ! call newG(S, 2, 1)
    ! call print_matrix(S%Gup, stdout, "newG computed G(2) = ")
    ! call print_matrix(A - S%Gup, stdout, "main - newG computed G(2) = ")






    contains


        subroutine make_A(S, A, l, sigma)
            type(Simulation), intent(in)  :: S
            real(dp)        , intent(out) :: A(S%N, S%N)
            integer         , intent(in)  :: l
            integer         , intent(in)  :: sigma

            integer :: i

            A = 0.0_dp
            do i = 1, S%N
                A(i, i) = exp(sigma * S%alpha * S%h(i, l))
            enddo


        endsubroutine make_A











endprogram main