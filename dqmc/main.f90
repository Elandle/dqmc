program main
    use numbertypes
    use simulate_mod
    implicit none
    
    type(Simulation) :: S
    integer          :: N           = 10
    integer          :: L           = 32
    integer          :: nstab       = 2
    integer          :: nbin        = 32
    integer          :: nmeassweep  = 32 * 10
    integer          :: nskip       = 2
    integer          :: nequil      = 50
    real(dp)         :: dtau        = 0.125_dp
    real(dp)         :: temp        = 10.0_dp
    real(dp)         :: U           = 4.0_dp
    real(dp)         :: mu          = 0.0_dp
    call setup_simulation(S, N, L, nstab, nbin, nmeassweep, nskip, nequil, &
                          dtau, temp, U, mu, "ckb.txt")
    S%h = 1
    call simulate(S)


endprogram main