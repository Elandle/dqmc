module measurements_mod
    use numbertypes
    use simulationsetup_mod
    use statistics_mod
    implicit none
    !
    ! Contains procedures for performing measurements during a simulation,
    ! and for doing statistics on those measurements.
    !
    ! All measurements are assumed to be binned, and measurement procedures should
    ! have the following form:
    !
    !       subroutine measure_quantity(S, i)
    !
    ! where i is the bin entry to fill.
    !
    ! For each quantity measured, S should contain a single bin workspace,
    ! which holds the data of a single bin before doing statistics on the bin.
    ! Additionally, for each quantity measured, S should contain a workspace for
    ! holding bin averages.
    !
    ! The main simulation subroutine (subroutine simulate, from simulate_mod)
    ! calls the subroutine measure in this module whenever a measurement is
    ! performed. Here should be a call to each measure_quantity subroutine as
    ! described above that is desired to be measured. Once all single bins are
    ! filled, it calls the avgbin subroutine. Here a method of filling a bin
    ! average should be included for each measurement performed. Once the simulation
    ! is done, a call to dostatistics is done, where the collected binned data
    ! should have statistics done on them to prepare them for printing collected
    ! results.
    !
    ! Note that the avgbin subroutine accounts for the effect of the sign problem
    ! on collecting measurements: the average of a bin should be divided by the
    ! average of the sign of determinants (Metropolis ratios) that was present
    ! when each measurement was being collected.
    !
    contains

        subroutine measure(S, i)
            !
            ! Called by the main simulation subroutine (subroutine simulate,
            ! from simulate_mod) whenever a measurement is to be performed.
            !
            ! Measurements are performed by calling a series of measurement
            ! procedures (which should be placed in this module), which
            ! are assumed to fill the ith bin workspace slot with a
            ! new measurement.
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            call measure_sgn(S, i)
            call measure_den(S, i)


        endsubroutine measure

    
        subroutine measure_den(S, i)
            !
            ! Fills the ith up and down density bin slots with the currently
            ! measured value of up and down densities.
            !
            ! This subroutine averages these densities over all sites in the lattice
            !
            ! Mathematically, the density on site j is:
            !
            ! den(j) = 1 - G(j, j)
            !
            ! where this density is evaluated for a specific spin (up or dn spin).
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            integer :: j

            S%updenbin(i) = 0
            S%dndenbin(i) = 0
            do j = 1, S%N
                S%updenbin(i) = S%updenbin(i) + 1.0_dp - S%Gup(j, j)
                S%dndenbin(i) = S%dndenbin(i) + 1.0_dp - S%Gdn(j, j)
            enddo
            S%updenbin(i) = S%sgn * S%updenbin(i) / S%N
            S%dndenbin(i) = S%sgn * S%dndenbin(i) / S%N


        endsubroutine measure_den


        subroutine measure_sgn(S, i)
            !
            ! Fills the ith sgn bin slot with the currently measured value of sgn.
            ! sgn being the sign of the product of determinants/Metropolis ratios:
            !
            ! sgn = sgn(Rup * Rdn)
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            S%sgnbin(i) = S%sgn


        endsubroutine measure_sgn


        subroutine avgbin(S, i)
            !
            ! Averages the ith bin measurement, storing the averages so the bin workspace
            ! can be used again, and accounts for the sign problem on averages by
            ! dividing these averages by the average of the sign of the determinant
            ! obtained when the bin was measured.
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            S%sgnbinavgs  (i) = real(sum(S%sgnbin), dp) / S%binsize
            S%updenbinavgs(i) = vector_avg(S%updenbin, S%binsize) / S%sgnbinavgs(i)
            S%dndenbinavgs(i) = vector_avg(S%dndenbin, S%binsize) / S%sgnbinavgs(i)

            
        endsubroutine avgbin


        subroutine dostatistics(S)
            !
            ! Computes the average and error of final obtained binned measurements
            ! using the jacknife method.
            !
            type(Simulation), intent(inout) :: S

            call jackknife(S%sgnbinavgs  , S%nbin, S%sgnavg  , S%sgnerr)
            call jackknife(S%updenbinavgs, S%nbin, S%updenavg, S%updenerr)
            call jackknife(S%dndenbinavgs, S%nbin, S%dndenavg, S%dndenerr)


        endsubroutine dostatistics


endmodule measurements_mod