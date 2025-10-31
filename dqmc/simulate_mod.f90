module simulate_mod
    use numbertypes
    use simulationsetup_mod
    use equalgreens_mod
    use metropolisratios_mod
    use measurements_mod
    use statistics_mod
    use iso_fortran_env, only: output_unit
    implicit none
    !
    ! Contains procedures for carrying out a DQMC simulation in full,
    ! and for completing sweeps.
    !
    contains

        !> \brief Performs one Monte Carlo sweep.
        !!
        !! Sweeps each imaginary time slice \f$l = 1, \dots, L\f$, performing
        !! one Monte Carlo sweep.
        subroutine sweep(S)
            type(Simulation), intent(inout) :: S

            integer :: l

            ! Sweep through imaginary time
            do l = 1, S%L

                ! Update Green's functions for this time slice
                call timeupdate(S, l, 1)
                call timeupdate(S, l, -1)

                ! Sweep through sites of the lattice at slice l
                call sweepslice(S, l)
                
            enddo
        endsubroutine sweep


        !> \brief Sweeps through a slice.
        !!
        !! Sweeps through all sites \f$i = 1, \dots, N\f$, proposing
        !! a flip to \f$h(i, l)\f$ at a fixed imaginary time slice \f$l\f$.
        subroutine sweepslice(S, l)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            integer  :: i
            real(dp) :: rand

            ! Sweep through all sites i = 1, ..., N at a fixed imaginary time step l
            do i = 1, S%N
                ! Propose flipping h(i, l)

                ! Calculate the acceptance probability R of that flip
                call greens_R(S, i, l, 1)
                call greens_R(S, i, l, -1)

                ! The probability of accepting the flip is R = |Rup * Rdn|
                ! The sign is needed for measurements
                S%upsgn = sgn(S%Rup)
                S%dnsgn = sgn(S%Rdn)
                S%sgn   = S%upsgn * S%dnsgn
                S%R     = S%sgn * S%Rup * S%Rdn

                ! Generate a random number uniformly between 0 and 1
                call random_number(rand)

                ! Metropolis algorithm
                if (rand .lt. S%R) then
                    ! Accept the flip

                    ! Flip h(i, l)
                    S%h(i, l) = -S%h(i, l)

                    ! Update the Green's functions
                    call flipupdate(S, i, 1)
                    call flipupdate(S, i, -1)

                else
                    ! Reject the flip
                endif
            enddo
        endsubroutine sweepslice


        !> \brief Performs warmup/equilibriation sweeps and
        !! gets things ready for calling (extra setup needed the very first sweep).
        subroutine warmup(S)
            type(Simulation), intent(inout) :: S

            integer :: i, l

            ! Very first sweep is slightly different since Green's functions
            ! have not been computed yet

            ! l = 1 imaginary time sweep:
            l = 1
            write(output_unit, "(a, i5, a)") "Warmup sweep ", 1, "..."
            write(S%dunit    , "(a, i5, a)") "Warmup sweep ", 1, "..."
            call newG(S, l, 1)
            call newG(S, l, -1)
            call sweepslice(S, l)

            do l = 2, S%L
                ! Sweep through imaginary time

                ! Update Green's functions for this time slice
                call timeupdate(S, l, 1)
                call timeupdate(S, l, -1)

                ! Sweep through sites of the lattice at slice l
                call sweepslice(S, l)
            enddo

            do i = 2, S%nequil
                write(output_unit, "(a, i5, a)") "Warmup sweep ", i, "..."
                write(S%dunit    , "(a, i5, a)") "Warmup sweep ", i, "..."
                call sweep(S)
            enddo
        endsubroutine warmup


        subroutine simulate(S)
            !
            ! Runs a DQMC simulation as specified by S
            !
            type(Simulation), intent(inout) :: S

            integer     :: i, j, k

            ! Print information about how many sweeps there are to the terminal
            write(output_unit, "(i5, a)") S%ntotal, " Total  sweeps"
            write(output_unit, "(i5, a)") S%nequil, " Warmup sweeps"
            write(output_unit, "(i5, a)") (S%nskip + 1) * S%binsize, " Sweeps per bin (approximately)"
            write(output_unit, "(i5, a)") S%nbin, " Bins"

            ! For debugging, prints information about the starting simulation
            ! call debug_setup_print(S)

            ! Warmup
            call warmup(S)

            ! First measurement sweep right after warmup
            write(output_unit, "(a, i5, a, i5, a)")     "Bin ", 1, " measured sweep    ", 1, " ..."
            call sweep(S)
            call measure(S, 1)
            

            ! Separate first bin loop since a measurement sweep has already been done
            do i = 2, S%binsize
                do j = 1, S%nskip
                    ! Nonmeasured sweeps
                    write(output_unit, "(a, i5, a, i5, a)") "Bin ", 1, " nonmeasured sweep ", j, " ..."
                    call sweep(S)
                enddo
                ! Measurement sweep
                write(output_unit, "(a, i5, a, i5, a)")     "Bin ", 1, " measured sweep    ", i, " ..."
                call sweep(S)
                call measure(S, i)
            enddo
            call avgbin(S, 1)


            ! Loop over bins
            do k = 2, S%nbin
                do i = 1, S%binsize
                    do j = 1, S%nskip
                        ! Nonmeasured sweeps
                        write(output_unit, "(a, i5, a, i5, a)") "Bin ", k, " nonmeasured sweep ", j, " ..."
                        call sweep(S)
                    enddo
                    ! Measurement sweep
                    write(output_unit, "(a, i5, a, i5, a)")     "Bin ", k, " measured sweep    ", i, " ..."
                    call sweep(S)
                    call measure(S, i)
                enddo
                call avgbin(S, k)
            enddo
            
            ! Do statistics
            call dostatistics(S)

            ! Output results
            call output(S)
        endsubroutine simulate


        subroutine debug_setup_print(S)
            type(Simulation) :: S

            integer :: i, j

            write(S%dunit, "(a, f17.8)") "Chemical potential mu    = ", S%mu
            write(S%dunit, "(a, f17.8)") "Time discritization dtau = ", S%dtau
            write(S%dunit, "(a, i5)")    "Imaginary time steps L   = ", S%L
            write(S%dunit, "(a, f17.8)") "beta = L * dtau          = ", S%beta
            write(S%dunit, "(a)") "Hopping matrix T = "
            call print_matrix(S%T, S%dunit)
            write(S%dunit, "(a)") "Hoppings:"
            write(S%dunit, "(a)") "i, j, T(i, j)"
            do i = 1, S%N
                do j = 1, S%N
                    if (abs(S%T(i, j)) .ge. 10e-4) then
                        write(S%dunit, "(i5, a, i5, f17.8)") i, " , ", j, S%T(i, j)
                    endif
                enddo
            enddo

            write(S%dunit    , "(i5, a)") S%ntotal, " Total  sweeps"
            write(S%dunit    , "(i5, a)") S%nequil, " Warmup sweeps"
            write(S%dunit    , "(i5, a)") (S%nskip + 1) * S%binsize, " Sweeps per bin (approximately)"
            write(S%dunit    , "(i5, a)") S%nbin, " Bins"

        endsubroutine debug_setup_print


        subroutine output(S)
            type(Simulation) :: S

            complex(dp) :: a, b, c

            integer :: i, j

            open(newunit=S%ounit, file=S%outfilename, action="write", status="replace")

            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average sign                               = ", S%sgnavg         , " +- ", S%sgnerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average total density                      = ", S%totaldenavg    , " +- ", S%totaldenerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average upden                              = ", S%updenavg       , " +- ", S%updenerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average dnden                              = ", S%dndenavg       , " +- ", S%dndenerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average KE                                 = ", S%kineticavg     , " +- ", S%kineticerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average CHEE                               = ", S%chemicalavg    , " +- ", S%chemicalerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average PE                                 = ", S%potentialavg   , " +- ", S%potentialerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average E                                  = ", S%energyavg      , " +- ", S%energyerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average antiferromagnetic structure factor = ", S%antiferroavg   , " +- ", S%antiferroerr
            write(S%ounit, "(a, 2f17.8, a, 2f17.8)") "Average uppol                              = ", S%uppolavg       , " +- ", S%uppolerr
            write(S%ounit, "(a, 2f17.8, a, 2f17.8)") "Average dnpol                              = ", S%dnpolavg       , " +- ", S%dnpolerr
            a = -(S%L**2)/(((2*4*atan(1.0_dp))**2)*S%N)
            b = a * log(S%uppolavg)**2
            c = abs(2 * a * log(b) / b) * S%uppolerr
            write(S%ounit, "(a, 2f17.8, a, 2f17.8)") "Average uplambda^2                         = ", b, " +- ", c
            b = a * log(S%dnpolavg)**2
            c = abs(2 * a * log(b) / b) * S%dnpolerr
            write(S%ounit, "(a, 2f17.8, a, 2f17.8)") "Average dnlambda^2                         = ", b, " +- ", c
            write(S%ounit, "(a)")                    "Gup                                        = "
            call print_matrix(S%Gupavg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_matrix(S%Guperr, S%ounit)
            write(S%ounit, "(a)")                    "Gdn                                        = "
            call print_matrix(S%Gdnavg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_matrix(S%Gdnerr, S%ounit)
            write(S%ounit, "(a)")                    "Spin density correlation                   = "
            call print_matrix(S%spindenscorravg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_matrix(S%spindenscorrerr, S%ounit)
            write(S%ounit, "(a)")                    "Spin spin    correlation                   = "
            call print_matrix(S%spinspincorravg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_matrix(S%spinspincorrerr, S%ounit)
            ! call abc(S)
            write(S%ounit, "(a)")                    "Average upden (full)                       = "
            call print_vector(S%updenfullavg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_vector(S%updenfullerr, S%ounit)
            write(S%ounit, "(a)")                    "Average dnden (full)                       = "
            call print_vector(S%dndenfullavg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_vector(S%dndenfullerr, S%ounit)
            write(S%ounit, "(a)")                    "Average double occupancy (full)            = "
            call print_vector(S%doubleoccfullavg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_vector(S%doubleoccfullerr, S%ounit)
            write(S%ounit, "(a)")                    "Average magnetic moment (full)             = "
            call print_vector(S%magmomentavg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_vector(S%magmomenterr, S%ounit)

            write(S%ounit, "(a)") "Bipartite matrix (entry i, j = 1 means sites i and j are on the same sublattice, -1 different) = "
            call print_integer_matrix(S%bipartsgn, S%ounit)


        endsubroutine output

        subroutine abc(S)
            type(Simulation) :: S

            integer :: i, j
            real(dp) :: sum

            sum = 0.0_dp

            do i = 1, S%N
                do j = 1, S%N
                    sum = sum + S%spindenscorravg(i, j) * ( (-1) ** (i+j))
                enddo
            enddo

            write(output_unit, "(f17.8)") sum / S%N


        endsubroutine abc


        integer function sgn(x)
            !
            ! Returns the sign of x as an integer.
            !
            !     1 if x >= 0
            !    -1 if x <  0
            !
            real(dp), intent(in) :: x

            if (x .ge. 0.0_dp) then
                sgn = 1
            else
                sgn = -1
            endif


        endfunction sgn


endmodule simulate_mod