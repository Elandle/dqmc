program main
    use numbertypes
    use simulate_mod
    implicit none
    
    type(Simulation) :: S
    integer          :: N           = 16
    integer          :: L           = 32
    integer          :: nstab       = 5
    integer          :: nbin        = 4
    integer          :: nmeassweep  = 4 * 200
    integer          :: nskip       = 5
    integer          :: nequil      = 200
    real(dp)         :: dtau        = 0.125_dp
    real(dp)         :: temp        = 10.0_dp
    real(dp)         :: U           = 2.0_dp
    real(dp)         :: mu          = 0.0_dp
    character(len=100)    :: filename = "squareckb.txt"
    call setup_simulation(S, N, L, nstab, nbin, nmeassweep, nskip, nequil, &
                          dtau, temp, U, mu, filename)
    S%h = 1
    call simulate(S)

    print *, "abc"
endprogram main