module simulationsetup_mod
    use numbertypes
    use readinputfile_mod
    use checkerboard_mod
    use convenientla_mod
    implicit none
    !
    ! Contains the Simulation datatype, usually denoted by S, and procedures
    ! for setting it up for a simulation.
    !
    ! A user should call appropriate setup procedures to configure S,
    ! so it may be called to perform a DQMC simulation, during which
    ! computations done during the simulation are stored in S.
    !
    ! Being the main datatype of a DQMC simulation, most procedures
    ! are coded to take S in as an argument. This is unfortunate for testing
    ! purposes (eg, constructing a Green's function for a given Hubbard-Stratonovich
    ! field), but can usually be worked around be calling a Simulation initialization
    ! procedure with arbitrary values for unnecessary variables or defining
    ! a Simulation datatype and filling in values as desired, and then calling the
    ! procedure. A list of dependencies on S for calling various procedures for
    ! testing intends to be developed.
    !
    type :: Simulation
        integer :: N          ! Number of sites
        integer :: L          ! Number of imaginary time slices
        
        integer :: nstab      ! Every nstab flips, the Green's functions will be computed from scratch
        integer :: north      ! How many single particle propagator's (B matrices) can be multiplied together accurately before doing a QRP
        integer :: nbin       ! How many bins to put measurements into (must divide nmeassweep)
        integer :: nmeassweep ! How many sweeps will have measurements performed
        integer :: nskip      ! How many sweeps to skip between measurements
        integer :: ntotal     ! Total number of sweeps performed (nequil + (nmeassweep - 1) * nskip + nmeassweep)
        integer :: nequil     ! How many sweeps to equilibriate
        integer :: binsize    ! How many measurements fit in each bin (nmeassweep / nbin)

        real(dp) :: dtau      ! Trotterization parameter
        real(dp) :: beta      ! Inverse temperature (beta = 1 / temp)
        real(dp) :: temp      ! Temperature
        real(dp) :: mu        ! Chemical potential
        real(dp) :: U         ! Electron interaction strength
        real(dp) :: alpha     ! Part of Hubbard Stratonovich transformation (alpha = arccosh(exp(dtau * U / 2)))
        real(dp) :: aldmu     ! Avoid multiplying by exp(dtau * mu) so often (aldmu = alpha * exp(dtau * mu))
        real(dp) :: aldmuinv  ! Avoid multiplying by 1 / aldmu so often (aldmuinv = 1 / aldmu)

        real(dp) :: Rup       ! up spin Metropolis weight
        real(dp) :: deltaup   ! Intermediate variable for changing Gup in response to a h(i, l) --> -h(i, l) flip
        real(dp) :: Rdn       ! dn spin Metropolis weight
        real(dp) :: deltadn   ! Intermediate variable for changing Gdn in response to a h(i, l) --> -h(i, l) flip

        real(dp) :: R         ! Metropolis weight (R = Rup * Rdn)
        real(dp) :: rand      ! Uniform random number between 0 and 1

        integer :: upstabi    ! Counter for determining when to compute Gup from scratch
        integer :: dnstabi    ! Counter for determining when to compute Gdn from scratch

        character(len=100)    :: ckbfilename ! File name of checkerboard file to read in
        integer               :: ckbiounit   ! Input/output unit for the ckbfilename file

        character(len=100)    :: outfilename
        integer               :: ounit

        character(len=100)    :: debfilename
        integer               :: dunit

        real(dp), allocatable :: T(:, :)     ! Hopping matrix (N x N)
        integer , allocatable :: h(:, :)     ! Hubbard Stratonovich field (N x L)
        real(dp), allocatable :: Gup(:, :)   ! up spin equal time Green's function (N X N)
        real(dp), allocatable :: Gdn(:, :)   ! dn spin equal time Green's function (N x N)
        integer               :: upsgn       ! sign of Rup
        integer               :: dnsgn       ! sign of Rdn
        integer               :: sgn         ! sign of R = Rup * Rdn

        integer               :: info        ! info argument for LAPACK
        
        ! Workspaces:
        real(dp), allocatable :: ckbwork(:)     ! Checkerboard work vector (N)
        real(dp), allocatable :: flipwork(:, :) ! Green's function fast flip update work array (N x 2; holds 2 vectors)
        
        ! Workspaces for carrying out the ASvQRD algorithm (stable chain multiplication)
        real(dp), allocatable :: qrdB   (:, :)    ! B matrix for ASvQRD (N x N)
        real(dp), allocatable :: qrdQ   (:, :)    ! Q matrix for ASvQRD (N x N)
        real(dp), allocatable :: qrdtau (:)       ! tau vector for ASvQRD (N)
        real(dp), allocatable :: qrdR   (:, :)    ! R matrix for ASvQRD (N x N)
        real(dp), allocatable :: qrdD   (:)       ! D vector for ASvQRD (N)
        real(dp), allocatable :: qrdF   (:)       ! F vector for ASvQRD (N)
        real(dp), allocatable :: qrdT   (:, :)    ! T matrix for ASvQRD (N x N)
        real(dp), allocatable :: qrdmatwork(:, :) ! work matrix for ASvQRD (N x N)
        real(dp), allocatable :: qrdwork(:)       ! work vector for ASvQRD (qrdlwork)
        integer , allocatable :: qrdP   (:)       ! P vector for ASvQRD (N)
        integer , allocatable :: qrdI   (:)       ! I vector for ASvQRD (N)
        integer               :: qrdlwork         ! length of qrdwork (at least 3*N+1)

        ! Workspaces for inverting general matrices
        real(dp), allocatable :: invwork(:)    ! work vector for inverting matrices (lwork)
        integer , allocatable :: invP   (:)    ! Permutation vector for inverting matrices (N)
        integer               :: invlwork      ! length of invwork (at least N)

        ! Measurement variables
        ! A typical measurement should have:
        !     a bin array of length binsize, for holding the measurement of a single bin
        !     a binavgs array of length nbin, for holding the average of a single bin
        !     an avg real number, for holding the final average value of the measurement
        !     an err real number, for holding the final error value of the measurement
        integer , allocatable :: sgnbin(:)       ! bin array for sign measurements (binsize)
        real(dp), allocatable :: sgnbinavgs(:)   ! binavgs array for sign measurements (nbin)
        real(dp)              :: sgnavg          ! avg real number for sign measurements
        real(dp)              :: sgnerr          ! err real number for sign measurements

        real(dp), allocatable :: updenbin(:)     ! bin array for up spin density measurements (binsize) 
        real(dp), allocatable :: updenbinavgs(:) ! binavgs array for up spin density measurements (nbin)
        real(dp)              :: updenavg        ! avg real number for up spin density measurements
        real(dp)              :: updenerr        ! err real number for up spin density measurements

        real(dp), allocatable :: dndenbin(:)     ! bin array for dn spin density measurements (binsize)
        real(dp), allocatable :: dndenbinavgs(:) ! binavgs array for dn spin density measurements (nbin)
        real(dp)              :: dndenavg        ! avg real number for dn spin density measurements
        real(dp)              :: dndenerr        ! err real number for dn spin density measurements

        real(dp), allocatable :: kineticbin(:)     ! (binsize)
        real(dp), allocatable :: kineticbinavgs(:) ! (nbin)
        real(dp)              :: kineticavg        !
        real(dp)              :: kineticerr        !

        real(dp), allocatable :: potentialbin(:)     ! (binsize)
        real(dp), allocatable :: potentialbinavgs(:) ! (nbin)
        real(dp)              :: potentialavg        !
        real(dp)              :: potentialerr        !

        real(dp), allocatable :: energybin(:)     ! (binsize)
        real(dp), allocatable :: energybinavgs(:) ! (nbin)
        real(dp)              :: energyavg        !
        real(dp)              :: energyerr        !

        real(dp), allocatable :: spindenscorrbin(:, :, :)     ! (N x N x binsize)
        real(dp), allocatable :: spindenscorrbinavgs(:, :, :) ! (N x N x nbin)
        real(dp), allocatable :: spindenscorravg(:, :)        ! (N x N)
        real(dp), allocatable :: spindenscorrerr(:, :)        ! (N x N)

        complex(dp), allocatable :: polP(:)                   ! (N)
        real   (dp), allocatable :: polB(:, :)                ! (N x N)
        complex(dp), allocatable :: polBZ(:, :)               ! (N x N)
        complex(dp)              :: polpmeas

        complex(dp), allocatable :: uppolbin(:)               ! (binsize)
        complex(dp), allocatable :: dnpolbin(:)               ! (binsize)
        complex(dp), allocatable :: uppolbinavgs(:)           ! (nbin)
        complex(dp), allocatable :: dnpolbinavgs(:)           ! (nbin)
        complex(dp)              :: uppolavg
        complex(dp)              :: dnpolavg
        complex(dp)              :: uppolerr
        complex(dp)              :: dnpolerr

        real(dp), allocatable :: updenfullbin(:, :)           ! (N x binsize)
        real(dp), allocatable :: dndenfullbin(:, :)           ! (N x binsize)
        real(dp), allocatable :: updenfullbinavgs(:, :)       ! (N x nbin)
        real(dp), allocatable :: dndenfullbinavgs(:, :)       ! (N x nbin)
        real(dp), allocatable :: updenfullavg(:)              ! (N)
        real(dp), allocatable :: dndenfullavg(:)              ! (N)
        real(dp), allocatable :: updenfullerr(:)              ! (N)
        real(dp), allocatable :: dndenfullerr(:)              ! (N)

        real(dp), allocatable :: doubleoccfullbin(:, :)       ! (N x binsize)
        real(dp), allocatable :: doubleoccfullbinavgs(:, :)   ! (N x nbin)
        real(dp), allocatable :: doubleoccfullavg(:)          ! (N)
        real(dp), allocatable :: doubleoccfullerr(:)          ! (N)


        ! Used only when using bmult_mod instead of bmultexact_mod.
        ! That is, when B matrices are multiplied by using the approximate matrix exponential
        ! by the checkerboard method instead of the exact matrix exponential
        type(checkerboard)    :: ckb         ! Checkerboard (for multiplication by exp(dtau * T))
        type(checkerboard)    :: ckbinv      ! Inverse checkerboard (for multiplication by exp(-dtau * T) = inv(exp(dtau * T))

        ! Used only when using bmultexact_mod instead of bmult_mod.
        ! That is, when B matrices are multiplied by using the exact matrix exponential instead
        ! of checkerboard method
        real(dp), allocatable :: expT(:, :)      ! exp(dtau * T) (N x N)
        real(dp), allocatable :: expTinv(:, :)   ! inv(exp(dtau * T)) (N x N)


        ! TODO:
        ! When finished, consolidate workspaces to a minimum
        ! Query for optimal workspace sizes


    endtype Simulation


    contains

    subroutine setup_simulation(S, N, L, nstab, north, nbin, nmeassweep, nskip, nequil, &
        dtau, U, mu, ckbfilename, outfilename, debfilename)
        !
        ! Main way of setting up a simulation datatype S for use in simulation.
        ! After calling setup_simulation, simulate(S) (from simulate_mod) should immediately
        ! be callable to perform a DQMC simulation.
        !
        type(Simulation)  , intent(inout) :: S
        integer           , intent(in)    :: N
        integer           , intent(in)    :: L
        integer           , intent(in)    :: nstab
        integer           , intent(in)    :: north
        integer           , intent(in)    :: nbin
        integer           , intent(in)    :: nmeassweep
        integer           , intent(in)    :: nskip
        integer           , intent(in)    :: nequil
        real(dp)          , intent(in)    :: dtau
        real(dp)          , intent(in)    :: U
        real(dp)          , intent(in)    :: mu
        character(len=*)  , intent(in)    :: ckbfilename
        character(len=*)  , intent(in)    :: outfilename
        character(len=*)  , intent(in)    :: debfilename


        S%N          = N
        S%L          = L

        S%nstab      = nstab
        S%north      = north
        S%nbin       = nbin
        S%nmeassweep = nmeassweep
        S%nskip      = nskip
        S%nequil     = nequil
        S%binsize    = nmeassweep / nbin
        S%ntotal     = nequil + (nmeassweep - 1) * nskip + nmeassweep

        S%dtau     = dtau
        S%beta     = dtau * L
        S%temp     = 1 / S%beta
        S%U        = U
        S%mu       = mu
        S%alpha    = acosh(exp(dtau * U / 2.0_dp))
        S%aldmu    = S%alpha * exp(dtau * mu)
        S%aldmuinv = 1 / S%aldmu

        S%ckbfilename = ckbfilename
        S%outfilename = outfilename
        S%debfilename = debfilename

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
        allocate(S%qrdmatwork(N, N))
        allocate(S%qrdwork(S%qrdlwork))
        allocate(S%qrdP(N))
        allocate(S%qrdI(N))
        
        S%invlwork = 5 * N
        allocate(S%invwork(S%invlwork))
        allocate(S%invP(N))

        allocate(S%h(N, L))

        allocate(S%Gup(N, N))
        allocate(S%Gdn(N, N))

        allocate(S%sgnbin(S%binsize))
        allocate(S%sgnbinavgs(nbin))

        allocate(S%updenbin(S%binsize))
        allocate(S%updenbinavgs(nbin))

        allocate(S%dndenbin(S%binsize))
        allocate(S%dndenbinavgs(nbin))

        allocate(S%kineticbin(S%binsize))
        allocate(S%kineticbinavgs(nbin))

        allocate(S%potentialbin(S%binsize))
        allocate(S%potentialbinavgs(nbin))

        allocate(S%energybin(S%binsize))
        allocate(S%energybinavgs(nbin))


        allocate(S%spindenscorrbin(N, N, S%binsize))
        allocate(S%spindenscorrbinavgs(N, N, nbin))
        allocate(S%spindenscorravg(N, N))
        allocate(S%spindenscorrerr(N, N))

        allocate(S%polP(N))
        allocate(S%polB(N, N))
        allocate(S%polBZ(N, N))
        call make_P(S%polP, L, N, int(sqrt(real(N, dp))))

        allocate(S%uppolbin(S%binsize))
        allocate(S%dnpolbin(S%binsize))
        allocate(S%uppolbinavgs(S%nbin))
        allocate(S%dnpolbinavgs(S%nbin))

        allocate(S%updenfullbin(N, S%binsize))
        allocate(S%dndenfullbin(N, S%binsize))
        allocate(S%updenfullbinavgs(N, S%nbin))
        allocate(S%dndenfullbinavgs(N, S%nbin))
        allocate(S%updenfullavg(N))
        allocate(S%dndenfullavg(N))
        allocate(S%updenfullerr(N))
        allocate(S%dndenfullerr(N))

        allocate(S%doubleoccfullbin(N, S%binsize))
        allocate(S%doubleoccfullbinavgs(N, S%nbin))
        allocate(S%doubleoccfullavg(N))
        allocate(S%doubleoccfullerr(N))








        call read_ckb(S%ckb   , ckbfilename, S%ckbiounit,  dtau)
        call read_ckb(S%ckbinv, ckbfilename, S%ckbiounit, -dtau)

        ! Makes the random number generator give the same results every time
        ! TODO: implement seed
        call random_init(.true., .true.)

        S%upstabi = 0
        S%dnstabi = 0

        allocate(S%T(N, N))
        allocate(S%expT(N, N))
        allocate(S%expTinv(N, N))

        call read_ckbT(S%T, N, ckbfilename, S%ckbiounit, 1.0_dp)
        call expm(S%T, S%expT, S%dtau)
        call expm(S%T, S%expTinv, -1.0_dp*S%dtau)




    endsubroutine setup_simulation


    subroutine make_P(P, L, N, x)
        !
        ! Sets:
        !
        !       P = exp(2*pi*i/L) * P
        !
        ! where the P on the right-hand side of assignment is the geometry-dependent
        ! polarization (diagonal) matrix (stored as a vector) P.
        !
        ! Currently hardcoded for a square lattice with x sites in the x direction
        ! (number of sites in the y direction is not needed)
        !
        !
        ! max = x / 2 (floored half)
        !
        ! For even x, site distances have the pattern (example for x=6):
        !
        ! 1 --- 2 --- 3 --- 4 --- 5 --- 6
        ! 0     1     2     3     2     1
        !
        !
        ! For odd x, site distances have the pattern (example for x=7):
        !
        ! 1 --- 2 --- 3 --- 4 --- 5 --- 6 --- 7
        ! 0     1     2     3     3     2     1
        !
        complex(dp), intent(out) :: P(N)
        integer    , intent(in)  :: L
        integer    , intent(in)  :: N
        integer    , intent(in)  :: x

        integer  :: i, j, inc, max, min
        logical  :: odd
        real(dp) :: pi

        if (mod(x, 2) .eq. 0) then
            odd = .false.
        else
            odd = .true.
        endif

        max = x / 2
        min = 0
        inc = -1
        j   = min
        i   = 1

        if (N .eq. 1) then
            j = 1
            P(1) = complex(real(j, dp), 0.0_dp)
        else
            do
                P(i) = complex(real(j, dp), 0.0_dp)
                if (odd .and. (j .eq. max)) then
                    i = i + 1
                    P(i) = P(i-1)
                endif
                if ((j .eq. max) .and. (inc .eq. 1)) then
                    inc = -1
                elseif ((j .eq. min) .and. (inc .eq. -1)) then
                    inc = 1
                endif
                j = j + inc
                if (i .eq. N) then
                    exit
                endif
                i = i + 1
            enddo
        endif

        pi = 4 * atan(1.0_dp)
        P  = exp(complex(0.0_dp, 2.0_dp * pi / L)) * P


    endsubroutine make_P


    subroutine setup_simulation_input(S, fname, funit, ounit)
        type(Simulation), intent(inout) :: S
        character(len=*), intent(in)    :: fname
        integer         , intent(out)   :: funit
        integer         , intent(in)    :: ounit

        type(parameter_values) :: param_values
        type(parameter_set)    :: param_set

        call readinputfile(fname, funit, ounit, param_values, param_set)
        call printparams(param_values, ounit)
        call setup_simulation(S,                        &
                              param_values%N,           &
                              param_values%L,           &
                              param_values%nstab,       &
                              param_values%north,       &
                              param_values%nbin,        &
                              param_values%nmeassweep,  &
                              param_values%nskip,       &
                              param_values%nequil,      &
                              param_values%dtau,        &
                              param_values%U,           &
                              param_values%mu,          &
                              param_values%ckbfilename, &
                              param_values%outfilename, &
                              param_values%debfilename)


    endsubroutine setup_simulation_input


endmodule simulationsetup_mod