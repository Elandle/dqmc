!> \brief Contains procedures for multiplying by and creating
!! single-electron propagators \f$B_\sigma\f$.
!! 
!! Mathematically, the single-electron propagator \f$B_\sigma(l)\f$ from time step \f$l-1\f$
!! to \f$l\f$ with spin \f$\sigma\f$ (\f$1\f$ for \f$\uparrow\f$, \f$-1\f$ for \f$\downarrow\f$) is:
!! \f[\begin{aligned}
!! B_\sigma(l) &= \exp(\Delta\tau\mu)\text{diag}(\exp(\sigma\alpha h(:, l)))\exp(\Delta\tau T) \\
!! &= \text{diag}(\exp(\sigma\alpha h(:, l) + \Delta\tau\mu))\exp(\Delta\tau T)
!! \end{aligned}\f]
!!
!! There are two options for multiplying by \f$\exp(\Delta\tau T)\f$.
!! The first (and recommended) way is by use of the checkerboard method.
!! The second is by exact multiplication (slow).
!! To keep the code simple (avoid having two modules to include or not include),
!! the procedures for both methods are present in this module.
!! The procedures for the method used should have the names:
!!
!! `right_Bmult left_Bmult right_Binvmult left_Binvmult make_B`
!!
!! When not in use, to keep the code organized, the checkerboard methods should have the names:
!!
!! `right_Bmult_ckb left_Bmult_ckb right_Binvmult_ckb left_Binvmult make_B_ckb`
!!
!! And the exact routines:
!!
!! `right_Bmult_exact left_Bmult_exact right_Binvmult_exact left_Binvmult_exact make_B_exact`
!!
!! Eventually this choice will be moved to a compilation flag, but for the time being
!! to keep the code simpler this method is being used now.
module bmult_mod
    use numbertypes
    use checkerboard_mod
    use simulationsetup_mod
    use customla_mod
    implicit none
    contains
    

        !> \brief Updates \f$A = AB_\sigma(l)\f$ using the checkerboard method.
        !!
        !! Recall:
        !! \f[\begin{aligned}
        !! B_\sigma(l) &= \exp(\Delta\tau\mu)\text{diag}(\exp(\sigma\alpha h(:, l)))\exp(\Delta\tau T) \\
        !! &= \text{diag}(\exp(\sigma\alpha h(:, l) + \Delta\tau\mu))\exp(\Delta\tau T)
        !! \end{aligned}\f]
        !! First \f$A\f$ is updated by \f$\text{diag}(\exp(\sigma\alpha h(:, l) + \Delta\tau\mu))\f$
        !! by column updates.
        !! Then \f$A\f$ is updated by \f$\exp(\Delta\tau T)\f$ by column updates with the checkerboard method.
        !!
        !! \see Todo: ensure combining \f$+\Delta\tau\mu\f$ into the `diag_mult` works as expected.
        subroutine right_Bmult(S, A, l, sigma)
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            ! Todo: check this works as expected
            ! call right_diagmult(A, exp(sigma * S%alpha * S%h(:, l) + S%dtau * S%mu), S%N)
            ! call right_ckbmult(S%ckb, A, S%N, S%ckbwork)

            call right_diagmult(A, exp(sigma * S%alpha * S%h(:, l)), S%N)
            call right_ckbmult(S%ckb, A, S%N, S%ckbwork)
            A = exp(S%dtau * S%mu) * A
        endsubroutine right_Bmult


        !> \brief Updates \f$A = B_\sigma(l)A\f$ using the checkerboard method.
        !!
        !! Recall:
        !! \f[\begin{aligned}
        !! B_\sigma(l) &= \exp(\Delta\tau\mu)\text{diag}(\exp(\sigma\alpha h(:, l)))\exp(\Delta\tau T) \\
        !! &= \text{diag}(\exp(\sigma\alpha h(:, l) + \Delta\tau\mu))\exp(\Delta\tau T)
        !! \end{aligned}\f]
        !! First \f$A\f$ is updated by \f$\exp(\Delta\tau T)\f$ by row updates with the checkerboard method.
        !! Then \f$A\f$ is updated by \f$\text{diag}(\exp(\sigma\alpha h(:, l) + \Delta\tau\mu))\f$ by row updates.
        !!
        !! \see Todo: ensure combining \f$+\Delta\tau\mu\f$ into the `diag_mult` works as expected.
        subroutine left_Bmult(S, A, l, sigma)
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            ! Todo: check this works as expected
            ! call left_ckbmult(S%ckb, A, S%N, S%ckbwork)
            ! call left_diagmult(A, exp(sigma * S%alpha * S%h(:, l) + S%dtau * S%mu), S%N)

            call left_ckbmult(S%ckb, A, S%N, S%ckbwork)
            call left_diagmult(A, exp(sigma * S%alpha * S%h(:, l)), S%N)
            A = exp(S%dtau * S%mu) * A
        endsubroutine left_Bmult


        subroutine right_Binvmult(S, A, l, sigma)
            !
            ! Updates:
            !
            ! A = A * inv(B(l))
            !
            ! Where B(l) is the lth single particle propagator from imaginary
            ! time (l-1)*dtau to l*dtau with spin sigma (1 for up, -1 for dn)
            !
            ! As a matrix:
            !
            ! B(l) = exp(dtau * mu) * diag(exp(sigma * alpha * h(:, l))) * exp(dtau * T)
            !
            ! inv(B(l)) = inv(exp(dtau * mu)) * inv(exp(dtau * T)) * inv(diag(exp(sigma * alpha * h(:, l))))
            !           = exp(-dtau * mu)     * exp(-dtau * T)     * diag(exp(-sigma * alpha * h(:, l)))
            !           =                       exp(-dtau * T)     * diag(exp(-sigma * alpha * h(:, l) - dtau * mu))
            !
            ! First A is updated by exp(-dtau * T) by row updates with the checkerboard method
            ! Then A is updated by diag(exp(-sigma * alpha * h(:, l) - dtau * mu)) by row updates
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            ! call right_ckbmult(S%ckbinv, A, S%N, S%ckbwork)
            ! call right_diagmult(A, exp(-sigma * S%alpha * S%h(:, l) - S%dtau * S%mu), S%N)

            call right_ckbmult(S%ckbinv, A, S%N, S%ckbwork)
            call right_diagmult(A, exp(-sigma * S%alpha * S%h(:, l)), S%N)
            A = exp(-S%dtau * S%mu) * A

        endsubroutine right_Binvmult


        subroutine left_Binvmult(S, A, l, sigma)
            !
            ! Updates:
            !
            ! A = inv(B(l)) * A
            !
            ! Where B(l) is the lth single particle propagator from imaginary
            ! time (l-1)*dtau to l*dtau with spin sigma (1 for up, -1 for dn)
            !
            ! As a matrix:
            !
            ! B(l) = exp(dtau * mu) * diag(exp(sigma * alpha * h(:, l))) * exp(dtau * T)
            !
            ! inv(B(l)) = inv(exp(dtau * mu)) * inv(exp(dtau * T)) * inv(diag(exp(sigma * alpha * h(:, l))))
            !           = exp(-dtau * mu)     * exp(-dtau * T)     * diag(exp(-sigma * alpha * h(:, l)))
            !           =                       exp(-dtau * T)     * diag(exp(-sigma * alpha * h(:, l) - dtau * mu))
            !
            ! First A is updated by exp(-dtau * T) by row updates with the checkerboard method
            ! Then A is updated by diag(exp(-sigma * alpha * h(:, l))) by row updates
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            ! call left_diagmult(A, exp(-sigma * S%alpha * S%h(:, l) - S%dtau * S%mu), S%N)
            ! call left_ckbmult(S%ckbinv, A, S%N, S%ckbwork)

            call left_diagmult(A, exp(-sigma * S%alpha * S%h(:, l)), S%N)
            call left_ckbmult(S%ckbinv, A, S%N, S%ckbwork)
            A = A / exp(S%dtau * S%mu)


        endsubroutine left_Binvmult


        subroutine make_B(S, A, l, sigma)
            !
            ! Sets:
            !
            ! A = B(l)
            !
            ! Where B(l) is the lth single particle propagator from imaginary
            ! time (l-1)*dtau to l*dtau with spin sigma (1 for up, -1 for dn)
            !
            ! A is set to B(l) by first setting A to the identity matrix
            ! then multiplying A on the right by B(l) (in column-major languages,
            ! such as Fortran, right multiplication by B(l) is faster than left).
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(out)   :: A(S%N, S%N)
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            ! A = id
            call dlaset('A', S%N, S%N, 0.0_dp, 1.0_dp, A, S%N)
            ! A = A * B(l)
            call right_Bmult(S, A, l, sigma)


        endsubroutine make_B


        ! ----------------------------------------------------------------------------

        ! 2. Exact


        subroutine right_Bmult_exact(S, A, l, sigma)
            !
            ! Updates:
            !
            ! A = A * B(l)
            !
            ! Where B(l) is the lth single particle propagator from imaginary
            ! time (l-1)*dtau to l*dtau with spin sigma (1 for up, -1 for dn)
            !
            ! As a matrix:
            !
            ! B(l) = exp(dtau * mu) * diag(exp(sigma * alpha * h(:, l))) * exp(dtau * T)
            !
            ! First A is updated by diag(exp(sigma * S%alpha * S%h(:, l))) by column updates
            ! Then A is updated by exp(dtau * T) by column updates with the checkerboard method
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            call right_diagmult(A, exp(sigma * S%alpha * S%h(:, l)), S%N)
            call right_matmul(A, S%expT, S%N, S%qrdB)
            ! TODO: combine with diagmult exp, once certain this is working
            A = exp(S%dtau * S%mu) * A


        endsubroutine right_Bmult_exact


        subroutine left_Bmult_exact(S, A, l, sigma)
            !
            ! Updates:
            !
            ! A = B(l) * A
            !
            ! Where B(l) is the lth single particle propagator from imaginary
            ! time (l-1)*dtau to l*dtau with spin sigma (1 for up, -1 for dn)
            !
            ! As a matrix:
            !
            ! B(l) = exp(dtau * mu) * diag(exp(sigma * alpha * h(:, l))) * exp(dtau * T)
            !
            ! First A is updated by exp(dtau * T) by row updates with the checkerboard method
            ! Then A is updated by diag(exp(sigma * alpha * h(:, l))) by row updates
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            call left_matmul(A, S%expT, S%N, S%qrdB)
            call left_diagmult(A, exp(sigma * S%alpha * S%h(:, l)), S%N)
            A = exp(S%dtau * S%mu) * A


        endsubroutine left_Bmult_exact


        subroutine right_Binvmult_exact(S, A, l, sigma)
            !
            ! Updates:
            !
            ! A = A * inv(B(l))
            !
            ! Where B(l) is the lth single particle propagator from imaginary
            ! time (l-1)*dtau to l*dtau with spin sigma (1 for up, -1 for dn)
            !
            ! As a matrix:
            !
            ! B(l) = exp(dtau * mu) * diag(exp(sigma * alpha * h(:, l))) * exp(dtau * T)
            !
            ! inv(B(l)) = inv(exp(dtau * mu)) * inv(exp(dtau * T)) * inv(diag(exp(sigma * alpha * h(:, l))))
            !           = exp(-dtau * mu)     * exp(-dtau * T)     * diag(exp(-sigma * alpha * h(:, l)))
            !
            !
            ! First A is updated by exp(-dtau * T) by row updates with the checkerboard method
            ! Then A is updated by diag(exp(-sigma * alpha * h(:, l))) by row updates
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            call right_matmul(A, S%expTinv, S%N, S%qrdB)
            call right_diagmult(A, exp(-sigma * S%alpha * S%h(:, l)), S%N)
            A = exp(-S%dtau * S%mu) * A


        endsubroutine right_Binvmult_exact


        subroutine left_Binvmult_exact(S, A, l, sigma)
            !
            ! Updates:
            !
            ! A = inv(B(l)) * A
            !
            ! Where B(l) is the lth single particle propagator from imaginary
            ! time (l-1)*dtau to l*dtau with spin sigma (1 for up, -1 for dn)
            !
            ! As a matrix:
            !
            ! B(l) = exp(dtau * mu) * diag(exp(sigma * alpha * h(:, l))) * exp(dtau * T)
            !
            ! inv(B(l)) = inv(exp(dtau * mu)) * inv(exp(dtau * T)) * inv(diag(exp(sigma * alpha * h(:, l))))
            !           = exp(-dtau * mu)     * exp(-dtau * T)     * diag(exp(-sigma * alpha * h(:, l)))
            !
            !
            ! First A is updated by exp(-dtau * T) by row updates with the checkerboard method
            ! Then A is updated by diag(exp(-sigma * alpha * h(:, l))) by row updates
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            call left_diagmult(A, exp(-sigma * S%alpha * S%h(:, l)), S%N)
            call left_matmul(A, S%expTinv, S%N, S%qrdB)
            A = A / exp(S%dtau * S%mu)


        endsubroutine left_Binvmult_exact


        subroutine make_B_exact(S, A, l, sigma)
            !
            ! Sets:
            !
            ! A = B(l)
            !
            ! Where B(l) is the lth single particle propagator from imaginary
            ! time (l-1)*dtau to l*dtau with spin sigma (1 for up, -1 for dn)
            !
            ! A is set to B(l) by first setting A to the identity matrix
            ! then multiplying A on the right by B(l) (in column-major languages,
            ! such as Fortran, right multiplication by B(l) is faster than left).
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(out)   :: A(S%N, S%N)
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            ! A = id
            call dlaset('A', S%N, S%N, 0.0_dp, 1.0_dp, A, S%N)
            ! A = A * B(l)
            call right_Bmult(S, A, l, sigma)


        endsubroutine make_B_exact















endmodule bmult_mod