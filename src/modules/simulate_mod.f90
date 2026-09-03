module simulate_mod
    use stduse
    use simulationsetup_mod
    use equalgreens_mod
    use metropolisratios_mod
    use measurements_mod
    use statistics_mod
    use utilities
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

            integer  :: l
            real(dp) :: updiff, dndiff
            logical  :: upcheck, dncheck

            ! Sweep through imaginary time
            do l = 1, S%L

                ! Update Green's functions for this time slice
                call timeupdate(S, l, 1, updiff, upcheck)
                call timeupdate(S, l, -1, dndiff, dncheck)

                if (S%stab%tuning) then
                    if (upcheck .and. dncheck) then
                        call stabupdate(S%stab, max(updiff, dndiff))
                    endif
                endif

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
            integer :: rupsgn, rdnsgn

            ! Sweep through all sites i = 1, ..., N at a fixed imaginary time step l
            do i = 1, S%N
                ! Propose flipping h(i, l)

                ! Calculate the acceptance probability R of that flip
                call greens_R(S, i, l, 1)
                call greens_R(S, i, l, -1)

                rupsgn = sgn(S%Rup)
                rdnsgn = sgn(S%Rdn)
                S%R    = abs(S%Rup * S%Rdn)

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

                    S%upsgn = S%upsgn * rupsgn
                    S%dnsgn = S%dnsgn * rdnsgn
                    S%sgn   = S%upsgn * S%dnsgn

                else
                    ! Reject the flip
                endif
            enddo
        endsubroutine sweepslice

        !> \brief Performs warmup/equilibriation sweeps and
        !! gets things ready for calling (extra setup needed the very first sweep).
        subroutine warmup(S)
            type(Simulation), intent(inout) :: S

            integer  :: i, l
            real(dp) :: updiff, dndiff
            logical  :: upcheck, dncheck

            ! Very first sweep is slightly different since Green's functions
            ! have not been computed yet

            ! l = 1 imaginary time sweep:
            l = 1
            write(stdout, "(a, i5, a)") "Warmup sweep ", 1, "..."
            call newG(S, l, 1)
            call newG(S, l, -1)
            call sweepslice(S, l)

            do l = 2, S%L
                ! Sweep through imaginary time

                ! Update Green's functions for this time slice
                call timeupdate(S, l, 1, updiff, upcheck)
                call timeupdate(S, l, -1, dndiff, dncheck)

                if (S%stab%tuning) then
                    if (upcheck .and. dncheck) then
                        call stabupdate(S%stab, max(updiff, dndiff))
                    endif
                endif

                ! Sweep through sites of the lattice at slice l
                call sweepslice(S, l)
            enddo

            do i = 2, S%nequil
                write(stdout, "(a, i5, a)") "Warmup sweep ", i, "..."
                call sweep(S)
            enddo

            
            write(stdout, "(a, i4)") "Final nstab = ", S%stab%n
            write(stdout, "(a, es14.4)") "Maximum dG seen = ", S%stab%maxdiff
            S%stab%tuning = .false.
        endsubroutine warmup

        subroutine simulate(S)
            !
            ! Runs a DQMC simulation as specified by S
            !
            type(Simulation), intent(inout) :: S

            integer     :: i, j, k

            ! Print information about how many sweeps there are to the terminal
            write(stdout, "(i5, a)") S%ntotal, " Total  sweeps"
            write(stdout, "(i5, a)") S%nequil, " Warmup sweeps"
            write(stdout, "(i5, a)") (S%nskip + 1) * S%binsize, " Sweeps per bin (approximately)"
            write(stdout, "(i5, a)") S%nbin, " Bins"

            ! For debugging, prints information about the starting simulation
            ! call debug_setup_print(S)

            ! Warmup
            call warmup(S)

            ! First measurement sweep right after warmup
            write(stdout, "(a, i5, a, i5, a)")     "Bin ", 1, " measured sweep    ", 1, " ..."
            call sweep(S)
            call measure(S, 1)
            

            ! Separate first bin loop since a measurement sweep has already been done
            do i = 2, S%binsize
                do j = 1, S%nskip
                    ! Nonmeasured sweeps
                    write(stdout, "(a, i5, a, i5, a)") "Bin ", 1, " nonmeasured sweep ", j, " ..."
                    call sweep(S)
                enddo
                ! Measurement sweep
                write(stdout, "(a, i5, a, i5, a)")     "Bin ", 1, " measured sweep    ", i, " ..."
                call sweep(S)
                call measure(S, i)
            enddo
            call avgbin(S, 1)


            ! Loop over bins
            do k = 2, S%nbin
                do i = 1, S%binsize
                    do j = 1, S%nskip
                        ! Nonmeasured sweeps
                        write(stdout, "(a, i5, a, i5, a)") "Bin ", k, " nonmeasured sweep ", j, " ..."
                        call sweep(S)
                    enddo
                    ! Measurement sweep
                    write(stdout, "(a, i5, a, i5, a)")     "Bin ", k, " measured sweep    ", i, " ..."
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
            integer     :: i, j

            open(newunit=S%ounit, file=S%outfilename, action="write", status="replace")
            write(S%ounit, "(a, a)")          "Summary       = ", S%summary
            write(S%ounit, "(a, i0)")         "seed          = ", S%seed
            write(S%ounit, "(a, i10)")        "N             = ", S%N
            write(S%ounit, "(a, i10)")        "L             = ", S%L
            write(S%ounit, "(a, f17.8)")      "dtau          = ", S%dtau
            write(S%ounit, "(a, f17.8)")      "beta          = ", S%L * S%dtau
            write(S%ounit, "(a, f17.8)")      "U             = ", S%U
            write(S%ounit, "(a, f17.8)")      "mu            = ", S%mu
            write(S%ounit, "(a, i10)")        "nstab         = ", S%nstab
            write(S%ounit, "(a, i10)")        "north         = ", S%north
            write(S%ounit, "(a, i10)")        "nbin          = ", S%nbin
            write(S%ounit, "(a, i10)")        "nmeassweep    = ", S%nmeassweep
            write(S%ounit, "(a, i10)")        "nskip         = ", S%nskip
            write(S%ounit, "(a, i10)")        "nequil        = ", S%nequil
            write(S%ounit, "(a, a)")          "ckbupfilename = ", S%ckbupfilename
            write(S%ounit, "(a, a)")          "ckbdnfilename = ", S%ckbdnfilename
            write(S%ounit, "(a, a)")          "bipfilename   = ", S%bipfilename
            write(S%ounit, "(a, a)")          "outfilename   = ", S%outfilename
            write(S%ounit, "(a, a)")          "debfilename   = ", S%debfilename
            write(S%ounit, "(a)")             "up hopping matrix i j Tup(i, j) Tup(j, i) ="
            do i = 1, S%N
                do j = i + 1, S%N
                    if (abs(S%Tup(i,j)) > 1.0e-12_dp .or. abs(S%Tup(j,i)) > 1.0e-12_dp) then
                        write(S%ounit, "(2i7, 2f17.8)") i, j, S%Tup(i,j), S%Tup(j,i)
                    endif
                enddo
            enddo
            write(S%ounit, "(a)")             "dn hopping matrix i j Tdn(i, j) Tdn(j, i) ="
            do i = 1, S%N
                do j = i + 1, S%N
                    if (abs(S%Tdn(i,j)) > 1.0e-12_dp .or. abs(S%Tdn(j,i)) > 1.0e-12_dp) then
                        write(S%ounit, "(2i7, 2f17.8)") i, j, S%Tdn(i,j), S%Tdn(j,i)
                    endif
                enddo
            enddo


            write(S%ounit, "(a)")             "Bipartite matrix S ="
            do i = 1, S%N
                do j = 1, S%N
                    write(S%ounit, "(i4)", advance="no") S%bipartsgn(i, j)
                enddo
                write(S%ounit, "(a)") ""
            enddo

            write(S%ounit, "(a)") ""

            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average sign                               = ", S%sgnavg         , " +- ", S%sgnerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average total density                      = ", S%totaldenavg    , " +- ", S%totaldenerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average upden                              = ", S%updenavg       , " +- ", S%updenerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average dnden                              = ", S%dndenavg       , " +- ", S%dndenerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average KE                                 = ", S%kineticavg     , " +- ", S%kineticerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average CHEE                               = ", S%chemicalavg    , " +- ", S%chemicalerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average PE                                 = ", S%potentialavg   , " +- ", S%potentialerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average E                                  = ", S%energyavg      , " +- ", S%energyerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average antiferromagnetic structure factor = ", S%antiferroavg   , " +- ", S%antiferroerr
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
            close(S%ounit)
        endsubroutine output

endmodule simulate_mod