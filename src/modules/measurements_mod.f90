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
    use stduse
    use simulationsetup_mod
    use statistics_mod
    use utilities
    implicit none

    type scamesr(nbin, binsize)
        integer, len :: nbin
        integer, len :: binsize

        real(dp) :: bin(binsize)
        real(dp) :: binavgs(nbin)
        real(dp) :: avg
        real(dp) :: err
    endtype scamesr

    type vecmesr(n, nbin, binsize)
        integer, len :: n
        integer, len :: nbin
        integer, len :: binsize

        real(dp) :: bin(n, binsize)
        real(dp) :: binavgs(n, nbin)
        real(dp) :: avg(n)
        real(dp) :: err(n)
    endtype vecmesr

    type matmesr(m, n, nbin, binsize)
        integer, len :: m
        integer, len :: n
        integer, len :: nbin
        integer, len :: binsize

        real(dp) :: bin(m, n, binsize)
        real(dp) :: binavgs(m, n, nbin)
        real(dp) :: avg(m, n)
        real(dp) :: err(m, n)
    endtype matmesr

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
            call measure_spinspincorr(S, i)
            call measure_greens(S, i)
            call measure_den_full(S, i)
            call measure_doubleocc_full(S, i)
            call measure_kinetic(S, i)
            call measure_potential(S, i)
            call measure_chemical(S, i)
            call measure_energy(S, i)
            call measure_magmoment(S, i)
            call measure_antiferro(S, i)
        endsubroutine measure

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

            associate(N => S%N, updenbin => S%updenbin, dndenbin => S%dndenbin, totaldenbin =>S%totaldenbin, sgn => S%sgn)
                updenbin   (i) = 0.0_dp
                dndenbin   (i) = 0.0_dp
                totaldenbin(i) = 0.0_dp
                do j = 1, N
                    updenbin   (i) = updenbin   (i) + ni(S, j,  1)
                    dndenbin   (i) = dndenbin   (i) + ni(S, j, -1)
                    totaldenbin(i) = totaldenbin(i) + ni(S, j,  1) + ni(S, j, -1)
                enddo
                updenbin   (i) = sgn * updenbin   (i) / N
                dndenbin   (i) = sgn * dndenbin   (i) / N
                totaldenbin(i) = sgn * totaldenbin(i) / N
            endassociate
        endsubroutine measure_den

        subroutine measure_kinetic(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i


            real(dp) :: kinetic
            integer  :: j, k

            ! real(dp) :: kinmuup, kinmudn
            ! integer  :: j

            ! CHECKED ONCE TO GIVE SAME RESULT AS DOUBLE FOR LOOP
            ! ENSURE CORRECTNESS
            ! kinetic = sum(S%T * transpose(S%Gup)) + sum(S%T * transpose(S%Gdn))

            associate(N => S%N, T => S%T, kineticbin => S%kineticbin, sgn => S%sgn)
                kinetic = 0.0_dp
                do j = 1, N
                    do k = 1, N
                        kinetic = kinetic - T(j, k) * cdicj(S, j, k,  1) &
                                          - T(j, k) * cdicj(S, j, k, -1)
                    enddo
                enddo

                kineticbin(i) = sgn * kinetic
            endassociate
        endsubroutine measure_kinetic

        subroutine measure_chemical(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            real(dp) :: chemical
            integer  :: j

            associate(N => S%N, mu => S%mu, chemicalbin => S%chemicalbin, sgn => S%sgn)
                chemical = 0.0_dp
                do j = 1, N
                    chemical = chemical - mu * ni(S, j,  1)    &
                                        - mu * ni(S, j, -1)
                enddo

                chemicalbin(i) = sgn * chemical
            endassociate
        endsubroutine measure_chemical

        subroutine measure_potential(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            real(dp) :: potential
            integer  :: j

            associate(N => S%N, U => S%U, potentialbin => S%potentialbin, sgn => S%sgn)
                potential = 0.0_dp
                do j = 1, N
                    ! Particle-hole symmetric form
                    ! potential = potential + U * (                                             &
                    !                                          ni(S, j, 1) * ni(S, j, -1)       &
                    !                             - 0.5_dp  * (ni(S, j, 1) + ni(S, j, -1))      &
                    !                             + 0.25_dp                                     &
                    !                             )

                    ! Non-particle-hole symmetric form
                    potential = potential + U * ni(S, j, 1) * ni(S, j, -1)          
                enddo
                ! potential = potential * S%U
                potentialbin(i) = sgn * potential
            endassociate
        endsubroutine measure_potential

        subroutine measure_energy(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            associate(energybin => S%energybin, kineticbin => S%kineticbin, potentialbin => S%potentialbin, chemicalbin => S%chemicalbin, sgn => S%sgn)
                ! energybin(i) = sgn * kineticbin(i) + sgn * potentialbin(i)
                energybin(i) = kineticbin  (i)      &
                             + potentialbin(i)      &
                             + chemicalbin (i)
            endassociate
        endsubroutine measure_energy

        subroutine measure_den_full(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            integer :: j

            associate(N => S%N, updenfullbin => S%updenfullbin, dndenfullbin => S%dndenfullbin, sgn => S%sgn)
                do j = 1, N
                    updenfullbin(j, i) = sgn * ni(S, j,  1)
                    dndenfullbin(j, i) = sgn * ni(S, j, -1)
                enddo
            endassociate
        endsubroutine measure_den_full

        subroutine measure_doubleocc_full(S, i)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            integer :: j

            do j = 1, S%N
                S%doubleoccfullbin(j, i) = S%sgn * ninj(S, j, 1, j, -1)
            enddo
        endsubroutine measure_doubleocc_full

        subroutine measure_spindenscorr(S, k)
            !
            ! <mx(i)mx(j)> = (del(i, j) - Gup(j, i)) * Gdn(i, j)
            !              + (del(i, j) - Gdn(j, i)) * Gup(i, j)
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

        subroutine measure_spinspincorr(S, k)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: k

            integer :: i, j

            do j = 1, S%N
                do i = 1, S%N
                    S%spinspincorrbin(i, j, k) = S%sgn *                       &
                                                (                              &
                                                    ninj(S, i,  1, j,  1)      &
                                                -   ninj(S, i,  1, j, -1)      &
                                                -   ninj(S, i, -1, j,  1)      &
                                                +   ninj(S, i, -1, j, -1)      &
                                                )
                enddo
            enddo
        endsubroutine measure_spinspincorr

        subroutine measure_antiferro(S, k)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: k

            integer :: i, j

            S%antiferrobin(k) = 0.0_dp

            do i = 1, S%N
                do j = 1, S%N
                    S%antiferrobin(k) =  S%antiferrobin(k)            &
                                       + S%bipartsgn(i, j)            &
                                       *                              &
                                       (                              &
                                           ninj(S, i,  1, j,  1)      &
                                       -   ninj(S, i,  1, j, -1)      &
                                       -   ninj(S, i, -1, j,  1)      &
                                       +   ninj(S, i, -1, j, -1)      &
                                    )
                enddo
            enddo
            S%antiferrobin(k) = S%antiferrobin(k) * S%sgn / S%N
        endsubroutine measure_antiferro

        subroutine measure_greens(S, k)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: k

            integer :: i, j

            do j = 1, S%N
                do i = 1, S%N
                    S%Gupbin(i, j, k) = S%sgn * S%Gup(i, j)
                    S%Gdnbin(i, j, k) = S%sgn * S%Gdn(i, j)
                enddo
            enddo
        endsubroutine measure_greens

        subroutine measure_magmoment(S, k)
            ! \langle\sigma_z^2\rangle = \langle(n_{i\uparrow}-n_{i\downarrow})^2\rangle
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: k

            integer :: i

            do i = 1, S%N
                S%magmomentbin(i, k) = S%sgn * 0.75_dp *               &
                                      (                                &
                                            ninj(S, i,  1, i,  1)      &
                                      -     ninj(S, i,  1, i, -1)      &
                                      -     ninj(S, i, -1, i,  1)      &
                                      +     ninj(S, i, -1, i, -1)      &
                                      )
            enddo
        endsubroutine measure_magmoment

        real(dp) function ninj(S, i, si, j, sj)
            !
            ! \langle n_{i\sigma}n_{j\sigma'} \rangle
            ! si = spin \sigma ; sj = spin \sigma' ; 1 or -1
            !
            type(Simulation), intent(in) :: S
            integer         , intent(in) :: i
            integer         , intent(in) :: si
            integer         , intent(in) :: j
            integer         , intent(in) :: sj

            if     (si .ne. sj) then
                if     (si .eq. 1) then
                    ninj = (1.0_dp - S%Gup(i, i)) * (1.0_dp - S%Gdn(j, j))
                elseif (si .eq. -1) then
                    ninj = (1.0_dp - S%Gdn(i, i)) * (1.0_dp - S%Gup(j, j))
                endif
            elseif (si .eq. sj) then
                if (si .eq. 1) then
                    ninj = (1.0_dp - S%Gup(i, i)) * (1.0_dp - S%Gup(j, j)) + (del(i, j) - S%Gup(j, i)) * S%Gup(i, j)
                elseif (si .eq. -1) then
                    ninj = (1.0_dp - S%Gdn(i, i)) * (1.0_dp - S%Gdn(j, j)) + (del(i, j) - S%Gdn(j, i)) * S%Gdn(i, j)
                endif
            endif
        endfunction ninj

        real(dp) function ni(S, i, spin)
            !
            ! Returns \langle n_{i\sigma} \rangle
            ! spin = \sigma ; 1 or -1
            !
            type(Simulation), intent(in) :: S
            integer         , intent(in) :: i
            integer         , intent(in) :: spin

            if (spin .eq. 1) then
                ni = 1.0_dp - S%Gup(i, i)
            else
                ni = 1.0_dp - S%Gdn(i, i)
            endif
        endfunction ni

        real(dp) function cicdj(S, i, j, spin)
            !
            ! Returns \langle c_{i\sigma}c_{j\sigma}^\dagger \rangle
            !
            type(Simulation), intent(in) :: S
            integer         , intent(in) :: i
            integer         , intent(in) :: j
            integer         , intent(in) :: spin

            if (spin .eq. 1) then
                cicdj = S%Gup(i, j)
            else
                cicdj = S%Gdn(i, j)
            endif
        endfunction cicdj

        real(dp) function cdicj(S, i, j, spin)
            !
            ! Returns \langle c_{i\sigma}^\dagger c_{j\sigma} \rangle
            !
            type(Simulation), intent(in) :: S
            integer         , intent(in) :: i
            integer         , intent(in) :: j
            integer         , intent(in) :: spin

            if (spin .eq. 1) then
                cdicj = del(i, j) - S%Gup(j, i)
            else
                cdicj = del(i, j) - S%Gdn(j, i)
            endif
        endfunction cdicj

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

            S%sgnbinavgs     (i) = real(sum(S%sgnbin), dp)  / S%binsize
            S%updenbinavgs   (i) =  vector_avg(S%updenbin   , S%binsize) / S%sgnbinavgs(i)
            S%dndenbinavgs   (i) =  vector_avg(S%dndenbin   , S%binsize) / S%sgnbinavgs(i)
            S%totaldenbinavgs(i) =  vector_avg(S%totaldenbin, S%binsize) / S%sgnbinavgs(i)

            do k = 1, S%N
                do j = 1, S%N
                    S%spindenscorrbinavgs(j, k, i) = vector_avg(S%spindenscorrbin(j, k, :), S%binsize) / S%sgnbinavgs(i)
                enddo
            enddo

            do k = 1, S%N
                do j = 1, S%N
                    S%spinspincorrbinavgs(j, k, i) = vector_avg(S%spinspincorrbin(j, k, :), S%binsize) / S%sgnbinavgs(i)
                enddo
            enddo

            do k = 1, S%N
                do j = 1, S%N
                    S%Gupbinavgs(j, k, i) = vector_avg(S%Gupbin(j, k, :), S%binsize) / S%sgnbinavgs(i)
                    S%Gdnbinavgs(j, k, i) = vector_avg(S%Gdnbin(j, k, :), S%binsize) / S%sgnbinavgs(i)
                enddo
            enddo

            do k = 1, S%N
                S%doubleoccfullbinavgs(k, i) = vector_avg(S%doubleoccfullbin(k, :), S%binsize) / S%sgnbinavgs(i)
                S%updenfullbinavgs    (k, i) = vector_avg(S%updenfullbin    (k, :), S%binsize) / S%sgnbinavgs(i)
                S%dndenfullbinavgs    (k, i) = vector_avg(S%dndenfullbin    (k, :), S%binsize) / S%sgnbinavgs(i)
                S%magmomentbinavgs    (k, i) = vector_avg(S%magmomentbin    (k, :), S%binsize) / S%sgnbinavgs(i)
            enddo

            S%kineticbinavgs     (i) = vector_avg(S%kineticbin  , S%binsize) / S%sgnbinavgs(i)
            S%potentialbinavgs   (i) = vector_avg(S%potentialbin, S%binsize) / S%sgnbinavgs(i)
            S%energybinavgs      (i) = vector_avg(S%energybin   , S%binsize) / S%sgnbinavgs(i)
            S%antiferrobinavgs   (i) = vector_avg(S%antiferrobin, S%binsize) / S%sgnbinavgs(i)
            S%chemicalbinavgs    (i) = vector_avg(S%chemicalbin , S%binsize) / S%sgnbinavgs(i)
        endsubroutine avgbin

        subroutine dostatistics(S)
            !
            ! Computes the average and error of final obtained binned measurements
            ! using the jacknife method.
            !
            type(Simulation), intent(inout) :: S

            integer :: i, j

            call  jackknife(S%sgnbinavgs     , S%nbin, S%sgnavg     , S%sgnerr)
            call  jackknife(S%updenbinavgs   , S%nbin, S%updenavg   , S%updenerr)
            call  jackknife(S%dndenbinavgs   , S%nbin, S%dndenavg   , S%dndenerr)
            call  jackknife(S%totaldenbinavgs, S%nbin, S%totaldenavg, S%totaldenerr)

            do i = 1, S%N
                do j = 1, S%N
                    call jackknife(S%spindenscorrbinavgs(i, j, :), S%nbin, S%spindenscorravg(i, j), S%spindenscorrerr(i, j))
                enddo
            enddo

            do i = 1, S%N
                do j = 1, S%N
                    call jackknife(S%spinspincorrbinavgs(i, j, :), S%nbin, S%spinspincorravg(i, j), S%spinspincorrerr(i, j))
                enddo
            enddo

            do i = 1, S%N
                do j = 1, S%N
                    call jackknife(S%Gupbinavgs(i, j, :), S%nbin, S%Gupavg(i, j), S%Guperr(i, j))
                    call jackknife(S%Gdnbinavgs(i, j, :), S%nbin, S%Gdnavg(i, j), S%Gdnerr(i, j))
                enddo
            enddo

            do i = 1, S%N
                call jackknife(S%doubleoccfullbinavgs(i, :), S%nbin, S%doubleoccfullavg(i), S%doubleoccfullerr(i))
                call jackknife(S%updenfullbinavgs    (i, :), S%nbin, S%updenfullavg    (i), S%updenfullerr    (i))
                call jackknife(S%dndenfullbinavgs    (i, :), S%nbin, S%dndenfullavg    (i), S%dndenfullerr    (i))
                call jackknife(S%magmomentbinavgs    (i, :), S%nbin, S%magmomentavg    (i), S%magmomenterr    (i))
            enddo

            call jackknife(S%kineticbinavgs  , S%nbin, S%kineticavg  , S%kineticerr)
            call jackknife(S%potentialbinavgs, S%nbin, S%potentialavg, S%potentialerr)
            call jackknife(S%energybinavgs   , S%nbin, S%energyavg   , S%energyerr)
            call jackknife(S%antiferrobinavgs, S%nbin, S%antiferroavg, S%antiferroerr)
            call jackknife(S%chemicalbinavgs , S%nbin, S%chemicalavg , S%chemicalerr)
        endsubroutine dostatistics

endmodule measurements_mod