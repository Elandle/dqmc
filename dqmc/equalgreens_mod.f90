module equalgreens_mod
    use numbertypes
    use simulationsetup_mod
    use asvqrd_mod
    use bmult_mod
    implicit none


    contains


        subroutine upflipupdate(S, i)
            !
            ! Updates:
            !
            ! G = G + x * outer(G(:, i) - ei, G(i, :))
            !
            ! Where:
            !
            ! x           = delta / (1 + (1 - G(i, i)) * delta)
            !             = delta / R
            ! delta       = exp(sigma * alpha * (hnew(i, l) - hold(i, l))) - 1
            !             = exp(2 * sigma * alpha * h(i, l)) - 1
            ! sigma       = up (1) or dn (-1) spin
            ! R           = spin sigma ratio of weights for Metropolis decision
            ! ei          = ith Cartesian unit vector
            ! outer(a, b) = outer product of the vectors a and b
            ! G           = up or dn equal time Green's function
            ! i           = site index a spin was flipped at
            ! l           = imaginary time index a spin was flipped at (not
            !               needed to be supplied for this subroutine)
            !
            ! This subroutine is hardcoded for the up spin Green's function
            !
            ! The Hubbard Stratonovich field h is assumed to be flipped
            ! at position (i, l) already. That is:
            !
            ! hnew(i, l) =  h(i, l)
            ! hold(i, l) = -h(i, l)
            !
            ! The ratio of weights R (Rup) and factor delta (deltaup)
            ! are assumed to be calculated beforehand and stored in S
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            real(dp) :: x

            ! TODO:
            ! Implement delayed update

            x = S%deltaup / S%Rup

            ! Vector G(:, i) - 1
            call dcopy(S%N, S%Gup(1, i), 1, S%flipwork(1, 1), 1)
            S%flipwork(i, 1) = S%flipwork(i, 1) - 1.0_dp

            ! Vector G(i, :)
            call dcopy(S%N, S%Gup(i, 1), S%N, S%flipwork(1, 2), 1)

            ! Outer product
            call dger(S%N, S%N, x, S%flipwork(1, 1), 1, S%flipwork(1, 2), 1, S%Gup, S%N)


            ! No BLAS:

            ! integer :: j, k
            ! Vector G(:, i) - 1
            ! S%flipwork(:, 1) = S%Gup(:, i)
            ! S%flipwork(i, 1) = S%flipwork(i, 1) - 1.0_dp

            ! Vector G(i, :)
            ! S%flipwork(:, 2) = S%Gup(i, :)

            ! Outer product
            ! do k = 1, S%N
            !     do j = 1, S%N
            !         S%Gup(j, k) = S%Gup(j, k) + x * S%flipwork(j, 1) * S%flipwork(k, 2)
            !     enddo
            ! enddo


        endsubroutine upflipupdate


        subroutine dnflipupdate(S, i)
            !
            ! Updates:
            !
            ! G = G + x * outer(G(:, i) - ei, G(i, :))
            !
            ! Where:
            !
            ! x           = delta / (1 + (1 - G(i, i)) * delta)
            !             = delta / R
            ! delta       = exp(sigma * alpha * (hnew(i, l) - hold(i, l))) - 1
            !             = exp(2 * sigma * alpha * h(i, l)) - 1
            ! sigma       = up (1) or dn (-1) spin
            ! R           = spin sigma ratio of weights for Metropolis decision
            ! ei          = ith Cartesian unit vector
            ! outer(a, b) = outer product of the vectors a and b
            ! G           = up or dn equal time Green's function
            ! i           = site index a spin was flipped at
            ! l           = imaginary time index a spin was flipped at (not
            !               needed to be supplied for this subroutine)
            !
            ! This subroutine is hardcoded for the dn spin Green's function
            !
            ! The Hubbard Stratonovich field h is assumed to be flipped
            ! at position (i, l) already. That is:
            !
            ! hnew(i, l) =  h(i, l)
            ! hold(i, l) = -h(i, l)
            !
            ! The ratio of weights R (Rdn) and factor delta (deltadn)
            ! are assumed to be calculated beforehand and stored in S
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i

            real(dp) :: x

            ! TODO:
            ! Implement delayed update

            x = S%deltadn / S%Rdn

            ! Vector G(:, i) - 1
            call dcopy(S%N, S%Gdn(1, i), 1, S%flipwork(1, 1), 1)
            S%flipwork(i, 1) = S%flipwork(i, 1) - 1.0_dp

            ! Vector G(i, :)
            call dcopy(S%N, S%Gdn(i, 1), S%N, S%flipwork(1, 2), 1)

            ! Outer product
            call dger(S%N, S%N, x, S%flipwork(1, 1), 1, S%flipwork(1, 2), 1, S%Gdn, S%N)


            ! No BLAS:

            ! integer :: j, k
            ! Vector G(:, i) - 1
            ! S%flipwork(:, 1) = S%Gdn(:, i)
            ! S%flipwork(i, 1) = S%flipwork(i, 1) - 1.0_dp

            ! Vector G(i, :)
            ! S%flipwork(:, 2) = S%Gdn(i, :)

            ! Outer product
            ! do k = 1, S%N
            !     do j = 1, S%N
            !         S%Gdn(j, k) = S%Gdn(j, k) + x * S%flipwork(j, 1) * S%flipwork(k, 2)
            !     enddo
            ! enddo


        endsubroutine dnflipupdate


        subroutine upwrap(S, l)
            !
            ! Updates:
            !
            ! G = B(l) * G * inv(B(l))
            !
            ! For use in wrapping the equal-time Green's function G in the following cases:
            !
            ! When a sweep through sites at an imaginary time step l - 1 < L is done, the Green's
            ! function G is wrapped/updated to the next imaginary step l by the update:
            !
            ! G = B(l) * G * inv(B(l))
            !
            ! When a sweep through sites at the final imaginary time step L is done, the Green's
            ! function G is wrapped/updated for the next sweep (which starts at imaginary time
            ! step 1) by the update:
            !
            ! G = B(1) * G * inv(B(1))
            !
            ! This subroutine is hardcoded for the up spin Green's function
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            call right_Binvmult(S, S%Gup, l, 1)
            call left_Bmult(S, S%Gup, l, 1)


        endsubroutine upwrap


        subroutine dnwrap(S, l)
            !
            ! Updates:
            !
            ! G = B(l) * G * inv(B(l))
            !
            ! For use in wrapping the equal-time Green's function G in the following cases:
            !
            ! When a sweep through sites at an imaginary time step l - 1 < L is done, the Green's
            ! function G is wrapped/updated to the next imaginary step l by the update:
            !
            ! G = B(l) * G * inv(B(l))
            !
            ! When a sweep through sites at the final imaginary time step L is done, the Green's
            ! function G is wrapped/updated for the next sweep (which starts at imaginary time
            ! step 1) by the update:
            !
            ! G = B(1) * G * inv(B(1))
            !
            ! This subroutine is hardcoded for the dn spin Green's function
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            call right_Binvmult(S, S%Gdn, l, 1)
            call left_Bmult(S, S%Gdn, l, 1)


        endsubroutine dnwrap


        subroutine uptimeupdate(S, l)
            !
            ! Updates Gup for a new time slice iteration in a sweep
            !
            ! Updates Gup in one of the two ways:
            !
            ! if this is the nstabth time Gup is being updated since
            ! a new computation from scratch, then Gup is computed from scratch
            !
            ! otherwise Gup is updated from the previous time iteration by wrapping
            ! (even if the last iteration was in a separate sweep call)
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            if (mod(S%upstabi+1, S%nstab) .eq. 0) then
                call newGup(S, l)
                S%upstabi = 0
            else
                call upwrap(S, l)
                S%upstabi = S%upstabi + 1
            endif


        endsubroutine uptimeupdate


        subroutine dntimeupdate(S, l)
            !
            ! Updates Gdn for a new time slice iteration in a sweep
            !
            ! Updates Gdn in one of the two ways:
            !
            ! if this is the nstabth time Gdn is being updated since
            ! a new computation from scratch, then Gdn is computed from scratch
            !
            ! otherwise Gdn is updated from the previous time iteration by wrapping
            ! (even if the last iteration was in a separate sweep call)
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            if (mod(S%dnstabi+1, S%nstab) .eq. 0) then
                call newGdn(S, l)
                S%dnstabi = 0
            else
                call dnwrap(S, l)
                S%dnstabi = S%dnstabi + 1
            endif


        endsubroutine dntimeupdate


        subroutine newGup(S, l)
            !
            ! Computes (spin omitted):
            !
            ! G(l) = inv(id + B(l) * ... * B(1) * B(L) * ... * B(l+1))
            !
            ! Which is the equal time Green's function at imaginary time step l.
            !
            ! This subroutine is hardcoded to compute the spin up equal time Green's function, and
            ! stores the result in S%Gup.
            !
            ! In other words, this subroutine sets:
            !
            ! S%Gup = Gup(l)
            !
            ! where Gup(l) is computed from scratch.
            !
            ! G(l) is computed over two steps:
            !
            !      T = B(l) * ... * B(1) * B(L) * ... * B(l+1)
            !      G = inv(id + T) done stabley
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            ! qrdT = B(l) * ... * B(1) * B(L) * ... * B(l+1)
            call newGupchain(S, l)

            ! D(L) = Db * Ds decomposition, see ASvQRD algorithm
            ! D = Db, F = Ds
            call DbDs(S%qrdD, S%qrdF, S%N)

            ! G = inv(id + qrdT)
            if (mod(S%L, 2) .eq. 0) then                                                                ! Q information stored in Q
                call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)             ! Q = full Q
                call right_diaginvmult(S%qrdQ, S%qrdD, S%N)                                                 ! Q = Q * inv(Db), no transpose
                call left_diagmult(S%qrdT, S%qrdF, S%N)                                                     ! T = Ds * T
                call add_trans(S%qrdT, S%qrdQ, S%N)                                                         ! T = T + trans(Q)
                call invert(S%qrdT, S%N, S%invP, S%invwork, S%invlwork, S%info)                             ! T = inv(T)
                call dgemm('n', 't', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdQ, S%N, 0.0_dp, S%Gup, S%N)   ! G = T * trans(Q)
            else                                                                                        ! Q information stored in B
                call dorgqr(S%N, S%N, S%N, S%qrdB, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)             ! B = full Q
                call right_diaginvmult(S%qrdB, S%qrdD, S%N)                                                 ! B = B * inv(Db), no transpose
                call left_diagmult(S%qrdT, S%qrdF, S%N)                                                     ! T = Ds * T
                call add_trans(S%qrdT, S%qrdB, S%N)                                                         ! T = T + trans(B)
                call invert(S%qrdT, S%N, S%invP, S%invwork, S%invlwork, S%info)                             ! T = inv(T)
                call dgemm('n', 't', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdB, S%N, 0.0_dp, S%Gup, S%N)   ! G = T * trans(B)
            endif


        endsubroutine newGup


        subroutine newGdn(S, l)
            !
            ! Computes (spin omitted):
            !
            ! G(l) = inv(id + B(l) * ... * B(1) * B(L) * ... * B(l+1))
            !
            ! Which is the equal time Green's function at imaginary time step l.
            !
            ! This subroutine is hardcoded to compute the spin dn equal time Green's function, and
            ! stores the result in S%Gdn.
            !
            ! In other words, this subroutine sets:
            !
            ! S%Gdn = Gdn(l)
            !
            ! where Gdn(l) is computed from scratch.
            !
            ! G(l) is computed over two steps:
            !
            !      T = B(l) * ... * B(1) * B(L) * ... * B(l+1)
            !      G = inv(id + T) done stabley
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            ! qrdT = B(l) * ... * B(1) * B(L) * ... * B(l+1)
            call newGdnchain(S, l)

            ! D(L) = Db * Ds decomposition, see ASvQRD algorithm
            ! D = Db, F = Ds
            call DbDs(S%qrdD, S%qrdF, S%N)

            ! G = inv(id + qrdT)
            if (mod(S%L, 2) .eq. 0) then                                                                ! Q information stored in Q
                call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)             ! Q = full Q
                call right_diaginvmult(S%qrdQ, S%qrdD, S%N)                                                 ! Q = Q * inv(Db), no transpose
                call left_diagmult(S%qrdT, S%qrdF, S%N)                                                     ! T = Ds * T
                call add_trans(S%qrdT, S%qrdQ, S%N)                                                         ! T = T + trans(Q)
                call invert(S%qrdT, S%N, S%invP, S%invwork, S%invlwork, S%info)                             ! T = inv(T)
                call dgemm('n', 't', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdQ, S%N, 0.0_dp, S%Gdn, S%N)   ! G = T * trans(Q)
            else                                                                                        ! Q information stored in B
                call dorgqr(S%N, S%N, S%N, S%qrdB, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)             ! B = full Q
                call right_diaginvmult(S%qrdB, S%qrdD, S%N)                                                 ! B = B * inv(Db), no transpose
                call left_diagmult(S%qrdT, S%qrdF, S%N)                                                     ! T = Ds * T
                call add_trans(S%qrdT, S%qrdB, S%N)                                                         ! T = T + trans(B)
                call invert(S%qrdT, S%N, S%invP, S%invwork, S%invlwork, S%info)                             ! T = inv(T)
                call dgemm('n', 't', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdB, S%N, 0.0_dp, S%Gdn, S%N)   ! G = T * trans(B)
            endif


        endsubroutine newGdn


        subroutine newGupchain(S, l)
            !
            ! Sets:
            !
            ! S%qrdT = Bup(l) * ... * Bup(1) * Bup(L) * ... * Bup(l+1)
            !
            ! where the multiplication is done stably using the ASvQRD algorithm.
            !
            ! For implementation details, see the subroutine chainiteration
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            integer :: j

            ! S%qrdB = B(l)
            call make_B(S, S%qrdB, l, 1)

            ! j = 1 iteration
            call firstchainiteration(S%qrdB, S%qrdtau, S%qrdR, S%qrdD, S%qrdT, S%qrdI, S%qrdP, S%qrdwork, S%qrdlwork, S%info, S%N)

            do j = 2, S%L
                if (mod(j, 2) .eq. 0) then
                    call make_B(S, S%qrdQ, getj(j, S%L, l), 1)
                    call chainiteration(S%qrdQ, S%qrdB, S%qrdtau, S%qrdR, S%qrdD, S%qrdT, S%qrdI, S%qrdP, S%qrdwork, S%qrdlwork, &
                                        S%info, S%N)
                else
                    call make_B(S, S%qrdB, getj(j, S%L, l), 1)
                    call chainiteration(S%qrdB, S%qrdQ, S%qrdtau, S%qrdR, S%qrdD, S%qrdT, S%qrdI, S%qrdP, S%qrdwork, S%qrdlwork, &
                                        S%info, S%N)
                endif
            enddo


        endsubroutine newGupchain


        subroutine newGdnchain(S, l)
            !
            ! Sets:
            !
            ! S%qrdT = Bdn(l) * ... * Bdn(1) * Bdn(L) * ... * Bdn(l+1)
            !
            ! where the multiplication is done stably using the ASvQRD algorithm.
            !
            ! For implementation details, see the subroutine chainiteration
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l

            integer :: j

            ! S%qrdB = B(l)
            call make_B(S, S%qrdB, l, -1)

            ! j = 1 iteration
            call firstchainiteration(S%qrdB, S%qrdtau, S%qrdR, S%qrdD, S%qrdT, S%qrdI, S%qrdP, S%qrdwork, S%qrdlwork, S%info, S%N)

            do j = 2, S%L
                if (mod(j, 2) .eq. 0) then
                    call make_B(S, S%qrdQ, getj(j, S%L, l), -1)
                    call chainiteration(S%qrdQ, S%qrdB, S%qrdtau, S%qrdR, S%qrdD, S%qrdT, S%qrdI, S%qrdP, S%qrdwork, S%qrdlwork, &
                                        S%info, S%N)
                else
                    call make_B(S, S%qrdB, getj(j, S%L, l), -1)
                    call chainiteration(S%qrdB, S%qrdQ, S%qrdtau, S%qrdR, S%qrdD, S%qrdT, S%qrdI, S%qrdP, S%qrdwork, S%qrdlwork, &
                                        S%info, S%N)
                endif
            enddo


        endsubroutine newGdnchain


        function getj(i, N, l) result(j)
            !
            ! Returns the index of the ith B matrix from the right in the multiplication chain
            ! required when computing a single particle Green's function G from scratch
            !
            ! More precisely, the lth imaginary time step Green's function G(l) (of an omitted spin)
            ! is (note that L is N in the argument of this subroutine, two L/l arguments are not allowed
            ! for a single subroutine):
            !
            ! G(l) = inv(id + B(l) * ... * B(1) * B(L) * ... * B(l+1))
            !
            ! Which requires computing the matrix product:
            !
            !         l                  L-l            l + L-l = L
            ! B(l) * ... * B(1) * B(L) * ... * B(l+1)
            !
            ! in a stable manner. This done by the ASvQRD algorithm (see the subroutine chainiteration).
            !
            ! The 1st B matrix from the right is B(l+1)
            ! 2nd is B(l+2)
            ! ...
            ! L-lth is B(l+L-l = L)
            ! L-l+1th is B(1 = l-L+(L-l+1))
            ! L-l+2th is B(2 = l-L+(L-l+2))
            ! ...
            ! L-l+l=Lth is B(l = l-L+(L))
            !
            ! The ASvQRD algorithm is done by iterating through an integer i = 1, ... L, (ie, through
            ! the 1st rightmost matrix to the Lth rightmost matrix). This function ensures that the
            ! correct B matrix is chosen for a given i when calculating the matrix product needed to
            ! compute G.
            !
            integer, intent(in)  :: i
            integer, intent(in)  :: N
            integer, intent(in)  :: l

            integer :: j

            if (i .le. N - l) then
                j = l + i
            else
                j = l - N + i
            endif

        endfunction getj


endmodule equalgreens_mod