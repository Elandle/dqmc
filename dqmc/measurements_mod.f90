    !> \brief Contains procedures for performing measurements during a simulation,
    !! and for doing statistics on those measurements.
    !!
    !! All measurements are assumed to be binned, and measurement procedures should
    !! have the following form:
    !!
    !!       subroutine measure_quantity(S, i)
    !!
    !! where i is the bin entry to fill.
    !!
    !! For each quantity measured, S should contain a single bin workspace,
    !! which holds the data of a single bin before doing statistics on the bin.
    !! Additionally, for each quantity measured, S should contain a workspace for
    !! holding bin averages.
    !!
    !! The main simulation subroutine (subroutine simulate, from simulate_mod)
    !! calls the subroutine measure in this module whenever a measurement is
    !! performed. Here should be a call to each measure_quantity subroutine as
    !! described above that is desired to be measured. Once all single bins are
    !! filled, it calls the avgbin subroutine. Here a method of filling a bin
    !! average should be included for each measurement performed. Once the simulation
    !! is done, a call to dostatistics is done, where the collected binned data
    !! should have statistics done on them to prepare them for printing collected
    !! results.
    !!
    !! Note that the avgbin subroutine accounts for the effect of the sign problem
    !! on collecting measurements: the average of a bin should be divided by the
    !! average of the sign of determinants (Metropolis ratios) that was present
    !! when each measurement was being collected.
module measurements_mod
    use numbertypes
    use simulationsetup_mod
    use polarization_mod
    use statistics_mod
    use iso_fortran_env, only: terminal => output_unit
    implicit none


    contains


        !> \brief Called by the main simulation subroutine (subroutine simulate,
        !! from simulate_mod) whenever a measurement is to be performed.
        !!
        !! Measurements are performed by calling a series of measurement
        !! procedures (which should be placed in this module), which
        !! are assumed to fill the ith bin workspace slot with a
        !! new measurement.
        subroutine measure(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            call measure_sgn(S, i)
            call measure_den(S, i)
            call measure_spindenscorr(S, i)
            call measure_pol(S, i)
            call measure_den_full(S, i)
            call measure_doubleocc_full(S, i)
            call measure_kinetic(S, i)
            call measure_potential(S, i)
            call measure_energy(S, i)
        endsubroutine measure


        subroutine measure_pol(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            call measure_P(S, S%polB, S%polBZ, S%polpmeas,  1)
            S%uppolbin(i) = S%sgn * S%polpmeas
            call measure_P(S, S%polB, S%polBZ, S%polpmeas, -1)
            S%dnpolbin(i) = S%sgn * S%polpmeas
        endsubroutine measure_pol


        !> \brief Fills the ith up and down density bin slots with the currently
        !! measured value of up and down densities.
        !!
        !! This subroutine averages these densities over all sites in the lattice
        !!
        !! Mathematically, the density on site j is:
        !!
        !! \f[\rho_\sigma(j) = 1 - G_\sigma(j, j)\f]
        !!
        !! where this density is evaluated for a specific spin (up or dn spin).
        subroutine measure_den(S, i)
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


        subroutine measure_kinetic(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i


            real(dp) :: kinetic
            real(dp) :: kinmuup, kinmudn
            integer  :: j

            kinetic = sum(S%T * transpose(S%Gup)) + sum(S%T * transpose(S%Gdn))
            kinmuup = 0.0_dp; kinmudn = 0.0_dp
            do j = 1, S%N
                kinmuup = kinmuup + 1.0_dp - S%Gup(j, j) !  + S%U / 2.0_dp
                kinmudn = kinmudn + 1.0_dp - S%Gdn(j, j) !  + S%U / 2.0_dp
            enddo
            kinetic = kinetic - S%mu * kinmuup - S%mu * kinmudn

            S%kineticbin(i) = S%sgn * kinetic
        endsubroutine measure_kinetic


        subroutine measure_potential(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            real(dp) :: potential
            integer  :: j

            potential = 0.0_dp
            do j = 1, S%N
                ! potential = potential + (S%Gup(j, j) - 0.5_dp) * (S%Gdn(j, j) - 0.5_dp) * S%U
                potential = potential + (1.0_dp - S%Gup(j, j)) * (1.0_dp - S%Gdn(j, j))
            enddo
            potential = potential * S%U
            S%potentialbin(i) = S%sgn * potential
        endsubroutine measure_potential


        subroutine measure_energy(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            S%energybin(i) = S%sgn * S%kineticbin(i) + S%sgn * S%potentialbin(i)
        endsubroutine measure_energy


        subroutine measure_den_full(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            integer :: j

            do j = 1, S%N
                S%updenfullbin(j, i) = S%Sgn * (1.0_dp - S%Gup(j, j))
                S%dndenfullbin(j, i) = S%Sgn * (1.0_dp - S%Gdn(j, j))
            enddo
        endsubroutine measure_den_full


        subroutine measure_doubleocc_full(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            integer :: j

            do j = 1, S%N
                S%doubleoccfullbin(j, i) = S%Sgn * (1.0_dp - S%Gup(j, j)) * (1.0_dp - S%Gdn(j, j))
            enddo
        endsubroutine measure_doubleocc_full

        subroutine measure_spindenscorr(S, k)
            !
            ! < mx(i)mx(j)> = (del(i, j) - Gup(j, i)) * Gdn(i, j)
            !               + (del(i, j) - Gdn(j, i)) * Gup(i, j)
            !
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: k

            integer :: i, j

            ! spindenscorrbin(i, j, k)
            ! measurement between sites i and j in bin k

            do j = 1, S%N
                do i = 1, S%N
                    S%spindenscorrbin(i, j, k) = S%sgn * ((del(i, j) - S%Gup(j, i)) * S%Gdn(i, j) &
                                                       +  (del(i, j) - S%Gdn(j, i)) * S%Gup(i, j))
                enddo
            enddo
        endsubroutine measure_spindenscorr


        integer function del(i, j)
            integer, intent(in) :: i
            integer, intent(in) :: j

            if (i .eq. j) then
                del = 1
            else
                del = -1
            endif

        endfunction del


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

            integer :: j, k

            S%sgnbinavgs  (i) = real(sum(S%sgnbin), dp) / S%binsize
            S%updenbinavgs(i) =  vector_avg(S%updenbin, S%binsize) / S%sgnbinavgs(i)
            S%dndenbinavgs(i) =  vector_avg(S%dndenbin, S%binsize) / S%sgnbinavgs(i)
            S%uppolbinavgs(i) = zvector_avg(S%uppolbin, S%binsize) / S%sgnbinavgs(i)
            S%dnpolbinavgs(i) = zvector_avg(S%dnpolbin, S%binsize) / S%sgnbinavgs(i)

            do k = 1, S%N
                do j = 1, S%N
                    S%spindenscorrbinavgs(j, k, i) = vector_avg(S%spindenscorrbin(j, k, :), S%binsize) / S%sgnbinavgs(i)
                enddo
            enddo

            do k = 1, S%N
                S%doubleoccfullbinavgs(k, i) = vector_avg(S%doubleoccfullbin(k, :), S%binsize) / S%sgnbinavgs(i)
                S%updenfullbinavgs    (k, i) = vector_avg(S%updenfullbin    (k, :), S%binsize) / S%sgnbinavgs(i)
                S%dndenfullbinavgs    (k, i) = vector_avg(S%dndenfullbin    (k, :), S%binsize) / S%sgnbinavgs(i)
            enddo

            S%kineticbinavgs  (i) = vector_avg(S%kineticbin  , S%binsize) / S%sgnbinavgs(i)
            S%potentialbinavgs(i) = vector_avg(S%potentialbin, S%binsize) / S%sgnbinavgs(i)
            S%energybinavgs   (i) = vector_avg(S%energybin   , S%binsize) / S%sgnbinavgs(i)
            

        endsubroutine avgbin


        subroutine dostatistics(S)
            !
            ! Computes the average and error of final obtained binned measurements
            ! using the jacknife method.
            !
            type(Simulation), intent(inout) :: S

            integer :: i, j

            call  jackknife(S%sgnbinavgs  , S%nbin, S%sgnavg  , S%sgnerr)
            call  jackknife(S%updenbinavgs, S%nbin, S%updenavg, S%updenerr)
            call  jackknife(S%dndenbinavgs, S%nbin, S%dndenavg, S%dndenerr)
            call zjackknife(S%uppolbinavgs, S%nbin, S%uppolavg, S%uppolerr)
            call zjackknife(S%dnpolbinavgs, S%nbin, S%dnpolavg, S%dnpolerr)

            do i = 1, S%N
                do j = 1, S%N
                    call jackknife(S%spindenscorrbinavgs(i, j, :), S%nbin, S%spindenscorravg(i, j), S%spindenscorrerr(i, j))
                enddo
            enddo

            do i = 1, S%N
                call jackknife(S%doubleoccfullbinavgs(i, :), S%nbin, S%doubleoccfullavg(i), S%doubleoccfullerr(i))
                call jackknife(S%updenfullbinavgs    (i, :), S%nbin, S%updenfullavg    (i), S%updenfullerr    (i))
                call jackknife(S%dndenfullbinavgs    (i, :), S%nbin, S%dndenfullavg    (i), S%dndenfullerr    (i))
            enddo

            call jackknife(S%kineticbinavgs  , S%nbin, S%kineticavg  , S%kineticerr)
            call jackknife(S%potentialbinavgs, S%nbin, S%potentialavg, S%potentialerr)
            call jackknife(S%energybinavgs   , S%nbin, S%energyavg   , S%energyerr)


        endsubroutine dostatistics


endmodule measurements_mod