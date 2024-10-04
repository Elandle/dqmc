module simulate_mod
    use numbertypes
    use simulationsetup_mod
    use equalgreens_mod
    use metropolisratios_mod
    implicit none


    contains

        subroutine sweep(S)
            type(Simulation), intent(inout) :: S

            integer :: l

            do l = 1, S%L
                ! Sweep through imaginary time

                ! Update Green's functions for this time slice
                call uptimeupdate(S, l)
                call dntimeupdate(S, l)

                ! Sweep through sites of the lattice at slice l
                call sweepslice(S, l)
                
            enddo

        endsubroutine sweep

        subroutine sweepslice(S, l)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            integer :: i

            ! Sweep through sites at a fixed imaginary time
            do i = 1, S%N
                ! Proposal: flip h(i, l)

                ! Calculating the acceptance probability R of that flip
                call greens_Rup(S, i, l)
                call greens_Rdn(S, i, l)
                S%R = S%Rup * S%Rdn

                ! Generating a random number uniformly between 0 and 1
                call random_number(S%rand)

                ! Metropolis algorithm
                if (S%rand .lt. S%R) then
                    ! Accept the flip

                    ! Flip h(i, l)
                    S%h(i, l) = -S%h(i, l)

                    ! Update Green's functions
                    call upflipupdate(S, i)
                    call dnflipupdate(S, i)

                else
                    ! Reject the flip
                endif

            enddo

        endsubroutine sweepslice


        subroutine warmup(S)
            type(Simulation), intent(inout) :: S

            integer :: i, l

            ! Very first sweep is slightly different since Green's functions
            ! have not been computed yet

            ! l = 1 imaginary time sweep:
            l = 1
            call newGup(S, l)
            call newGdn(S, l)
            call sweepslice(S, l)

            do l = 2, S%L
                ! Sweep through imaginary time

                ! Update Green's functions for this time slice
                call uptimeupdate(S, l)
                call dntimeupdate(S, l)

                ! Sweep through sites of the lattice at slice l
                call sweepslice(S, l)
                
            enddo

            do i = 2, S%nequil
                call sweep(S)
            enddo

            
        endsubroutine warmup


        subroutine simulate(S)
            type(Simulation), intent(inout) :: S

            integer :: i, j, k


            ! Warmup
            call warmup(S)


            ! First measurement sweep right after warmup
            call sweep(S)
            ! call measure
            

            ! Separate first bin loop since a measurement sweep has already been done
            do i = 2, S%binsize
                do j = 1, S%nskip
                    ! Nonmeasured sweeps
                    call sweep(S)
                enddo
                ! Measurement sweep
                call sweep(S)
                ! call measure
            enddo


            ! Loop over bins
            do k = 2, S%nbin
                do i = 1, S%binsize
                    do j = 1, S%nskip
                        ! Nonmeasured sweeps
                        call sweep(S)
                    enddo
                    ! Measurement sweep
                    call sweep(S)
                    ! call measure
                enddo
            enddo
            
            
            ! Do statistics
            ! call statistics


            ! Output results
            ! call output


        endsubroutine simulate






endmodule simulate_mod