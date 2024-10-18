module simulate_mod
    use numbertypes
    use simulationsetup_mod
    use equalgreens_mod
    use metropolisratios_mod
    use measurements_mod
    use statistics_mod
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

            integer :: i

            ! Sweep through sites at a fixed imaginary time
            do i = 1, S%N
                ! Proposal: flip h(i, l)

                ! Calculating the acceptance probability R of that flip
                call greens_R(S, i, l, 1)
                call greens_R(S, i, l, -1)
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
                    call flipupdate(S, i, 1)
                    call flipupdate(S, i, -1)

                else
                    ! Reject the flip
                endif
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
            print *, "Warmup sweep ", 1, "..."
            call newG(S, l, 1)
            call newG(S, l, -1)
            call sweepslice(S, l)

            do l = 2, S%L
                ! Sweep through imaginary time

                ! Update Green's functions for this time slice
                ! call timeupdate(S, l, 1)
                ! call timeupdate(S, l, -1)
                call newG(S, l, 1)
                call newG(S, l, -1)

                ! Sweep through sites of the lattice at slice l
                call sweepslice(S, l)
                
            enddo

            do i = 2, S%nequil
                print *, "Warmup sweep ", i, "..."
                call sweep(S)
            enddo

            
        endsubroutine warmup


        subroutine simulate(S)
            !
            ! Runs a DQMC simulation as specified by S
            !
            type(Simulation), intent(inout) :: S

            integer :: i, j, k

            print *, S%ntotal, "total sweeps"
            print *, S%nequil, "warmup sweeps"
            print *, (S%nskip + 1) * S%binsize, "sweeps per bin (approximately)"
            print *, S%nbin, "bins"
            ! Warmup
            call warmup(S)


            ! First measurement sweep right after warmup
            print *, "Bin ", 1, "measured sweep ", 1, "..."
            call sweep(S)
            call measure(S, 1)
            

            ! Separate first bin loop since a measurement sweep has already been done
            do i = 2, S%binsize
                do j = 1, S%nskip
                    ! Nonmeasured sweeps
                    print *, "Bin ", 1, "nonmeasured sweep ", j, "..."
                    call sweep(S)
                enddo
                ! Measurement sweep
                print *, "Bin ", 1, "measured sweep...", i, "..."
                call sweep(S)
                call measure(S, i)
            enddo
            call avgbin(S, 1)


            ! Loop over bins
            do k = 2, S%nbin
                do i = 1, S%binsize
                    do j = 1, S%nskip
                        ! Nonmeasured sweeps
                        print *, "Bin ", k, "nonmeasured sweep ", j, "..."
                        call sweep(S)
                    enddo
                    ! Measurement sweep
                    print *, "Bin ", k, "measured sweep ", i, "..."
                    call sweep(S)
                    call measure(S, i)
                enddo
                call avgbin(S, k)
            enddo
            
            
            ! Do statistics
            call dostatistics(S)


            ! Output results
            ! call output
            print *, "Average sign  = ", S%sgnavg, "+-", S%sgnerr
            print *, "Average upden = ", S%updenavg, "+-", S%updenerr
            print *, "Average dnden = ", S%dndenavg, "+-", S%dndenerr


        endsubroutine simulate


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