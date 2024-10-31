module bmult_mod
    use numbertypes
    use checkerboard_mod
    use simulationsetup_mod
    use customla_mod
    implicit none
    !
    ! Contains procedures for multiplying by and creating
    ! single-electron propagators making use of the checkerboard method.
    !
    ! Mathematically, the single-electron propagator B(l) from time step l-1
    ! to l with spin sigma (1 for up, -1 for dn) is:
    !
    ! B(l) = exp(dtau * mu) * diag(exp(sigma * alpha * h(:, l)))             * exp(dtau * T)
    !      =                  diag(exp(sigma * alpha * h(:, l) + dtau * mu)) * exp(dtau * T)
    !
    ! Multiplication by exp(dtau * T) is done by use of the checkerboard method.
    !
    contains


        subroutine right_Bmult(S, A, l, sigma)
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
            ! B(l) = exp(dtau * mu) * diag(exp(sigma * alpha * h(:, l)))             * exp(dtau * T)
            !      =                  diag(exp(sigma * alpha * h(:, l) + dtau * mu)) * exp(dtau * T)
            !
            ! First A is updated by diag(exp(sigma * S%alpha * S%h(:, l) + dtau * mu)) by column updates
            ! Then A is updated by exp(dtau * T) by column updates with the checkerboard method
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            call right_diagmult(A, exp(sigma * S%alpha * S%h(:, l) + S%dtau * S%mu), S%N)
            call right_ckbmult(S%ckb, A, S%N, S%ckbwork)


        endsubroutine right_Bmult


        subroutine left_Bmult(S, A, l, sigma)
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
            ! B(l) = exp(dtau * mu) * diag(exp(sigma * alpha * h(:, l)))             * exp(dtau * T)
            !      =                  diag(exp(sigma * alpha * h(:, l) + dtau * mu)) * exp(dtau * T)
            !
            ! First A is updated by exp(dtau * T) by row updates with the checkerboard method
            ! Then A is updated by diag(exp(sigma * alpha * h(:, l) + dtau * mu)) * exp(dtau * T) by row updates
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(inout) :: A(S%N, S%N)
            integer         , intent(in)    :: sigma
            integer         , intent(in)    :: l

            call left_ckbmult(S%ckb, A, S%N, S%ckbwork)
            call left_diagmult(A, exp(sigma * S%alpha * S%h(:, l) + S%dtau * S%mu), S%N)


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

            call right_ckbmult(S%ckbinv, A, S%N, S%ckbwork)
            call right_diagmult(A, exp(-sigma * S%alpha * S%h(:, l) - S%dtau * S%mu), S%N)

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

            call left_diagmult(A, exp(-sigma * S%alpha * S%h(:, l) - S%dtau * S%mu), S%N)
            call left_ckbmult(S%ckbinv, A, S%N, S%ckbwork)


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


endmodule bmult_mod