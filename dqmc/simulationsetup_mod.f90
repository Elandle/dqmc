module simulationsetup_mod
    use numbertypes
    use checkerboard_mod
    implicit none


    type :: Simulation
        integer :: N          ! Number of sites
        integer :: L          ! Number of imaginary time slices
        
        integer :: nstab      ! Every nstab flips, the Green's functions will be computed from scratch
        integer :: nbin       ! How many bins to put measurements into (must divide nmeassweep)
        integer :: nmeassweep ! How many sweeps will have measurements performed
        integer :: nskip      ! How many sweeps to skip between measurements
        integer :: ntotal     ! Total number of sweeps performed (nequil + (nmeassweep - 1) * nskip + nmeassweep)
        integer :: nequil     ! How many sweeps to equilibriate
        integer :: binsize    ! How many measurements fit in each bin (nmeassweep / nbin)

        real(dp) :: dtau     ! Trotterization parameter
        real(dp) :: beta     ! Inverse temperature (beta = 1 / temp)
        real(dp) :: temp     ! Temperature
        real(dp) :: mu       ! Chemical potential
        real(dp) :: U        ! Electron interaction strength
        real(dp) :: alpha    ! Part of Hubbard Stratonovich transformation (arccosh(exp(dtau * U / 2)))
        real(dp) :: aldmu    ! Avoid multiplying by exp(dtau * mu) so often (alpha * exp(dtau * mu))
        real(dp) :: aldmuinv ! Avoid multiplying by 1 / aldmu so often (1 / aldmu)

        real(dp) :: Rup
        real(dp) :: deltaup
        real(dp) :: Rdn
        real(dp) :: deltadn

        real(dp) :: R
        real(dp) :: rand

        character(len=100) :: ckbfilename

        integer :: upstabi ! Counter for determining when to compute Gup from scratch
        integer :: dnstabi

        real(dp), allocatable :: T(:, :)   ! Hopping matrix (N x N)
        type(checkerboard)    :: ckb       ! Checkerboard (for multiplication by exp(dtau * T))
        type(checkerboard)    :: ckbinv    ! Inverse checkerboard (for multiplication by exp(-dtau * T) = inv(exp(dtau * T))
        integer , allocatable :: h(:, :)   ! Hubbard Stratonovich field (N x L)
        real(dp), allocatable :: Gup(:, :) ! Up spin equal time Green's function (N X N)
        real(dp), allocatable :: Gdn(:, :) ! Dn spin equal time Green's funciton (N x N)

        integer               :: info      ! info argument for LAPACK
        
        ! Workspaces:
        real(dp), allocatable :: ckbwork(:)     ! Checkerboard work vector (N)
        real(dp), allocatable :: flipwork(:, :) ! Green's function fast flip update work array (N x 2; holds 2 vectors)
        
        ! Workspaces for carrying out the ASvQRD algorithm (stable chain multiplication)
        real(dp), allocatable :: qrdB   (:, :) ! B matrix for ASvQRD (N x N)
        real(dp), allocatable :: qrdQ   (:, :) ! Q matrix for ASvQRD (N x N)
        real(dp), allocatable :: qrdtau (:)    ! tau vector for ASvQRD (N)
        real(dp), allocatable :: qrdR   (:, :) ! R matrix for ASvQRD (N x N)
        real(dp), allocatable :: qrdD   (:)    ! D vector for ASvQRD (N)
        real(dp), allocatable :: qrdF   (:)    ! F vector for ASvQRD (N)
        real(dp), allocatable :: qrdT   (:, :)    ! T matrix for ASvQRD (N x N)
        real(dp), allocatable :: qrdwork(:)    ! work vector for ASvQRD (qrdlwork)
        integer , allocatable :: qrdP   (:)    ! P vector for ASvQRD (N)
        integer , allocatable :: qrdI   (:)    ! I vector for ASvQRD (N)
        integer               :: qrdlwork      ! length of qrdwork (at least 3*N+1)

        ! Workspaces for inverting general matrices
        real(dp), allocatable :: invwork(:)    ! work vector for inverting matrices (lwork)
        integer , allocatable :: invP   (:)    ! P vector for inverting matrices (N)
        integer               :: invlwork      ! length of invwork (at least N)


        ! TODO:
        ! When finished, consolidate workspaces to a minimum
        ! Query for optimal workspace sizes


    endtype Simulation


    contains

    subroutine setup_simulation(S, N, L, nstab, nbin, nmeassweep, nskip, nequil, &
        dtau, temp, U, mu, ckbfilename)
        type(Simulation)  , intent(inout) :: S
        integer           , intent(in)    :: N
        integer           , intent(in)    :: L
        integer           , intent(in)    :: nstab
        integer           , intent(in)    :: nbin
        integer           , intent(in)    :: nmeassweep
        integer           , intent(in)    :: nskip
        integer           , intent(in)    :: nequil
        real(dp)          , intent(in)    :: dtau
        real(dp)          , intent(in)    :: temp
        real(dp)          , intent(in)    :: U
        real(dp)          , intent(in)    :: mu
        character(len=100), intent(in)    :: ckbfilename


        S%N          = N
        S%L          = L

        S%nstab      = nstab
        S%nbin       = nbin
        S%nmeassweep = nmeassweep
        S%nskip      = nskip
        S%nequil     = nequil
        S%binsize    = nmeassweep / nbin
        S%ntotal     = nequil + (nmeassweep - 1) * nskip + nmeassweep

        S%dtau     = dtau
        S%temp     = temp
        S%U        = U
        S%mu       = mu
        S%beta     = 1 / temp
        S%alpha    = acosh(exp(dtau * U / 2))
        S%aldmu    = S%alpha * exp(dtau * mu)
        S%aldmuinv = 1 / S%aldmu

        S%ckbfilename = ckbfilename

        S%qrdlwork = 5 * (3 * N + 1)
        allocate(S%ckbwork(N))
        allocate(S%flipwork(N, 2))
        allocate(S%qrdB(N, N))
        allocate(S%qrdQ(N, N))
        allocate(S%qrdtau(N))
        allocate(S%qrdR(N, N))
        allocate(S%qrdD(N))
        allocate(S%qrdF(N))
        allocate(S%qrdT(N, N))
        allocate(S%qrdwork(S%qrdlwork))
        allocate(S%qrdP(N))
        allocate(S%qrdI(N))
        
        S%invlwork = 5 * N
        allocate(S%invwork(S%invlwork))
        allocate(S%invP(N))

        allocate(S%h(N, L))


    endsubroutine setup_simulation


endmodule simulationsetup_mod