program main
    use numbertypes
    use simulate_mod
    implicit none
    
    type(Simulation) :: S
    integer          :: N           = 100
    integer          :: L           = 100
    integer          :: nstab       = 5
    integer          :: nbin        = 5
    integer          :: nmeassweep  = 5 * 1
    integer          :: nskip       = 1
    integer          :: nequil      = 5
    real(dp)         :: dtau        = 0.125_dp
    real(dp)         :: temp        = 10.0_dp
    real(dp)         :: U           = 4.0_dp
    real(dp)         :: mu          = 0.0_dp
    character(len=100)    :: filename = "squareckb.txt"
    call setup_simulation(S, N, L, nstab, nbin, nmeassweep, nskip, nequil, &
                          dtau, temp, U, mu, filename)
    S%h = 1
    call simulate(S)

    print *, "abc"
endprogram main