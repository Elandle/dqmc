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


        subroutine sweep(S)
            !
            ! Sweeps through all imaginary time slices 1, ..., L, performing
            ! one Monte Carlo sweep.
            !
            type(Simulation), intent(inout) :: S

            integer :: l

            do l = 1, S%L
                ! Sweep through imaginary time

                ! Update Green's functions for this time slice
                call timeupdate(S, l, 1)
                call timeupdate(S, l, -1)

                ! Sweep through sites of the lattice at slice l
                call sweepslice(S, l)
                
            enddo


        endsubroutine sweep

        subroutine sweepslice(S, l)
            !
            ! Sweeps slice l in a DQMC simulation.
            !
            ! Sweeps through all sites 1, ..., N, proposing a flip,
            ! for a given imaginary time slice l
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            ! Present for debugging ---
            real(dp) :: wr(S%N), wi(S%N)
            ! -------------------------

            integer :: i

            ! Sweep through sites at a fixed imaginary time
            do i = 1, S%N
                ! Proposal: flip h(i, l)

                ! Calculating the acceptance probability R of that flip
                ! Use one of the following:
                ! 1. From the Green's function (fast)
                ! call greens_R(S, i, l, 1)
                ! call greens_R(S, i, l, -1)
                ! 2. From the definition (slow)
                call scratch_R(S, i, l, 1)
                call scratch_R(S, i, l, -1)

                S%upsgn = sgn(S%Rup)
                S%dnsgn = sgn(S%Rdn)
                S%sgn   = S%upsgn * S%dnsgn
                S%R     = S%sgn * S%Rup * S%Rdn

                ! Generating a random number uniformly between 0 and 1
                call random_number(S%rand)

                ! Metropolis algorithm
                if (S%rand .lt. S%R) then
                    ! Accept the flip

                    ! Flip h(i, l)
                    S%h(i, l) = -S%h(i, l)

                    ! Update Green's functions
                    ! Use one of the following:
                    ! 1. Rank one update (fast)
                    ! call flipupdate(S, i, 1)
                    ! call flipupdate(S, i, -1)
                    ! 2. Make from scratch (slow)
                    call newG(S, l,  1)
                    call newG(S, l, -1)

                else
                    ! Reject the flip
                endif

                ! Present for debugging ---------
                write(S%dunit, "(a, i5, i5)") "Slice l, site i = ", l, i
                write(S%dunit, "(a)")     "Gup = "
                call print_matrix(S%Gup, S%dunit)
                call eigenvalues(S%Gup, wr, wi)
                write(S%dunit, "(a)") "Gup wr = "
                call print_vector(wr, S%dunit)
                write(S%dunit, "(a)") "Gup wi = "
                call print_vector(wi, S%dunit)

                write(S%dunit, "(a)")     "Gdn = "
                call print_matrix(S%Gdn, S%dunit)
                call eigenvalues(S%Gdn, wr, wi)
                write(S%dunit, "(a)") "Gdn wr = "
                call print_vector(wr, S%dunit)
                write(S%dunit, "(a)") "Gdn wi = "
                call print_vector(wi, S%dunit)
                ! -------------------------------

            enddo


        endsubroutine sweepslice


        subroutine warmup(S)
            !
            ! Performs the warmup sweeps of a DQMC simulation as specified by S
            !
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
            complex(dp) :: a, b, c

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

            write(output_unit, "(i5, a)") S%ntotal, " Total  sweeps"
            write(output_unit, "(i5, a)") S%nequil, " Warmup sweeps"
            write(output_unit, "(i5, a)") (S%nskip + 1) * S%binsize, " Sweeps per bin (approximately)"
            write(output_unit, "(i5, a)") S%nbin, " Bins"
            write(S%dunit    , "(i5, a)") S%ntotal, " Total  sweeps"
            write(S%dunit    , "(i5, a)") S%nequil, " Warmup sweeps"
            write(S%dunit    , "(i5, a)") (S%nskip + 1) * S%binsize, " Sweeps per bin (approximately)"
            write(S%dunit    , "(i5, a)") S%nbin, " Bins"

            ! Warmup
            call warmup(S)


            ! First measurement sweep right after warmup
            write(output_unit, "(a, i5, a, i5, a)")     "Bin ", 1, " measured sweep    ", 1, " ..."
            write(S%dunit    , "(a, i5, a, i5, a)")     "Bin ", 1, " measured sweep    ", 1, " ..."
            call sweep(S)
            call measure(S, 1)
            

            ! Separate first bin loop since a measurement sweep has already been done
            do i = 2, S%binsize
                do j = 1, S%nskip
                    ! Nonmeasured sweeps
                    write(output_unit, "(a, i5, a, i5, a)") "Bin ", 1, " nonmeasured sweep ", j, " ..."
                    write(S%dunit    , "(a, i5, a, i5, a)") "Bin ", 1, " nonmeasured sweep ", j, " ..."
                    call sweep(S)
                enddo
                ! Measurement sweep
                write(output_unit, "(a, i5, a, i5, a)")     "Bin ", 1, " measured sweep    ", i, " ..."
                write(S%dunit    , "(a, i5, a, i5, a)")     "Bin ", 1, " measured sweep    ", i, " ..."
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
                        write(S%dunit    , "(a, i5, a, i5, a)") "Bin ", k, " nonmeasured sweep ", j, " ..."
                        call sweep(S)
                    enddo
                    ! Measurement sweep
                    write(output_unit, "(a, i5, a, i5, a)")     "Bin ", k, " measured sweep    ", i, " ..."
                    write(S%dunit    , "(a, i5, a, i5, a)")     "Bin ", k, " measured sweep    ", i, " ..."
                    call sweep(S)
                    call measure(S, i)
                enddo
                call avgbin(S, k)
            enddo
            
            
            ! Do statistics
            call dostatistics(S)


            ! Output results
            ! call output

            open(newunit=S%ounit, file=S%outfilename, action="write", status="replace")

            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average sign                    = ", S%sgnavg      , " +- ", S%sgnerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average upden                   = ", S%updenavg    , " +- ", S%updenerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average dnden                   = ", S%dndenavg    , " +- ", S%dndenerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average KE                      = ", S%kineticavg  , " +- ", S%kineticerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average PE                      = ", S%potentialavg, " +- ", S%potentialerr
            write(S%ounit, "(a, f17.8, a, f17.8)")   "Average E                       = ", S%energyavg   , " +- ", S%energyerr
            write(S%ounit, "(a, 2f17.8, a, 2f17.8)") "Average uppol                   = ", S%uppolavg    , " +- ", S%uppolerr
            write(S%ounit, "(a, 2f17.8, a, 2f17.8)") "Average dnpol                   = ", S%dnpolavg    , " +- ", S%dnpolerr
            a = -(S%L**2)/(((2*4*atan(1.0_dp))**2)*S%N)
            b = a * log(S%uppolavg)**2
            c = abs(2 * a * log(b) / b) * S%uppolerr
            write(S%ounit, "(a, 2f17.8, a, 2f17.8)") "Average uplambda^2              = ", b, " +- ", c
            b = a * log(S%dnpolavg)**2
            c = abs(2 * a * log(b) / b) * S%dnpolerr
            write(S%ounit, "(a, 2f17.8, a, 2f17.8)") "Average dnlambda^2              = ", b, " +- ", c
            write(S%ounit, "(a)")                    "Spin density correlation        = "
            call print_matrix(S%spindenscorravg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_matrix(S%spindenscorrerr, S%ounit)
            ! call abc(S)
            write(S%ounit, "(a)")                    "Average upden (full)            = "
            call print_vector(S%updenfullavg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_vector(S%updenfullerr, S%ounit)
            write(S%ounit, "(a)")                    "Average dnden (full)            = "
            call print_vector(S%dndenfullavg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_vector(S%dndenfullerr, S%ounit)
            write(S%ounit, "(a)")                    "Average double occupancy (full) = "
            call print_vector(S%doubleoccfullavg, S%ounit)
            write(S%ounit, "(a)") "+-"
            call print_vector(S%doubleoccfullerr, S%ounit)
            

        endsubroutine simulate

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

        subroutine print_matrix(A, ounit)
            real(dp), intent(in) :: A(:, :)
            integer , intent(in) :: ounit
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            do i = 1, n
                do j = 1, m
                    write(ounit, "(f17.8)", advance="no") A(i, j)
                enddo
                write(ounit, "(a)") ""
            enddo


        endsubroutine print_matrix


        subroutine print_vector(v, ounit)
            real(dp), intent(in) :: v(:)
            integer , intent(in) :: ounit

            integer :: i, m

            m = size(v, 1)

            do i = 1, m
                write(ounit, "(f17.8)", advance="no") v(i)
            enddo
            write(ounit, "(a)") ""


        endsubroutine print_vector


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