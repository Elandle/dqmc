program main
    use numbertypes
    use simulate_mod
    implicit none
    !
    ! Executable file that runs a simulation.
    ! At the moment this is being used as a simulation tester,
    ! so it is heavily under construction. The eventual plan is to
    ! move over to a bunch of similar example files that run
    ! simulations like this one.
    !
    type(Simulation)   :: S
    integer            :: N           = 16
    integer            :: L           = 20
    integer            :: nstab       = 2
    integer            :: nbin        = 5
    integer            :: nmeassweep  = 5 * 100
    integer            :: nskip       = 5
    integer            :: nequil      = 100
    real(dp)           :: dtau        = 0.125_dp
    real(dp)           :: U           = 2.0_dp
    real(dp)           :: mu          = 0.0_dp
    character(len=100) :: filename = "squareckb.txt"



    
    call setup_simulation(S, N, L, nstab, nbin, nmeassweep, nskip, nequil, &
                          dtau, U, mu, filename)
    S%h = 1
    call simulate(S)






    contains


        subroutine print_matrix(A)
            real(dp), intent(in) :: A(:, :)
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            do j = 1, n
                do i = 1, m
                    write(*, "(F12.6)", advance="no") A(i, j)
                enddo
                write(*, *) ""
            enddo


        endsubroutine print_matrix


endprogram main