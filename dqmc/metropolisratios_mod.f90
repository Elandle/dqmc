module metropolisratios_mod
    use numbertypes
    use simulationsetup_mod
    use equalgreens_mod
    implicit none
    !
    ! Contains procedures for calculating the Metropolis ratio:
    !
    !     R = Rup * Rdn
    !
    ! that is used for determining whether or not a flip is accepted or not.
    !
    !
    ! TODO:
    ! Create subroutines to calculate R according to the formula:
    !
    ! R = det(id + B(beta, 0))
    !   = det(id + B(L) * ... * B(1))
    !
    ! using the ASvQRD algorithm.
    !
    ! This subroutine is not meant to be used in the final code, but rather
    ! testing purposes (ie, make sure R is calculated correctly)
    !
    contains


        subroutine greens_R(S, i, l, sigma)
            !
            ! Sets:
            !
            ! S%R     = 1 + (1 - G(i, i)) * delta
            ! S%delta = exp(-2 * sigma * alpha * h(i, l)) - 1
            !
            ! Where S%R is either S%Rup or S%Rdn and S%delta is either
            ! S%deltaup or S%deltadn, depending on whether sigma = 1 (spin up)
            ! or sigma = -1 (spin dn).
            !
            ! Note that this subroutine assumes that h(i, l) has not been flipped (since
            ! the purpose of this subroutine is to decide whether or not to flip h(i, l))
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            if (sigma .eq. 1) then      ! up spin (sigma = 1)
                S%deltaup = exp(-2 * sigma * S%alpha * S%h(i, l)) - 1
                S%Rup     = 1 + (1 - S%Gup(i, i)) * S%deltaup
            else                        ! dn spin (sigma = -1)
                S%deltadn = exp(-2 * sigma * S%alpha * S%h(i, l)) - 1
                S%Rdn     = 1 + (1 - S%Gdn(i, i)) * S%deltadn
            endif


        endsubroutine greens_R


        subroutine scratch_R(S, i, l, sigma)
            !
            ! Computes the Metropolis ratio for spin sigma (1 up, -1 dn)
            ! from the most basic definition in DQMC.
            !
            ! In more detail, the weight of a configuration of a Hubbard Stratonovich
            ! field h is:
            !
            !       w(h) = det(id + B(beta, 0))
            !            = det(id + B(L) * B(L-1) * ... * B(1))
            !
            ! A new h is proposed by h(i, l) = -h(i, l)
            ! (note: this code temporarily swaps S%h(i, l) = -S%h(i, l) but
            ! swaps it back at the end)
            !
            ! The Metropolis ratio is:
            !
            !       R = w(h') / w(h)
            !
            ! where h' is the proposed state.
            !
            ! This ratio is calculated by computing both w(h') and w(h)
            ! from determinants and then taking the ratio.
            !
            ! Note: the multiplication B(L) * ... * B(1) is done stably,
            ! but then the id is added and a determinant taken.
            ! Possible instability adding id.
            !
            ! For use in other routines, both S%R and S%delta are set.
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            integer :: j    ! B matrix counter
            integer :: ii   ! north counter
            integer :: info

            real(dp) :: Rold, Rnew


            ! Iteration j = 1
            j = 1

            ! Q = B(L)
            call make_B(S, S%qrdQ, getj(j, S%L, S%L), sigma)
                    
            if (j .lt. S%north) then
                do j = j+1, S%north
                    call left_Bmult(S, S%qrdQ, getj(j, S%L, S%L), sigma)
                enddo
            endif
            j = S%north

            ! QRP factorise Q
            S%qrdP = 0
            call dgeqp3(S%N, S%N, S%qrdQ, S%N, S%qrdP, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

            ! D = diag(Q)
            call diag(S%qrdQ, S%qrdD, S%N)

            ! T = uppertri(Q)
            call uppertri(S%qrdQ, S%qrdT, S%N)
            ! T = inv(D) * T
            call left_diaginvmult(S%qrdT, S%qrdD, S%N)
            ! T = T * inv(P)
            S%qrdI = S%qrdP
            call invert_permutation(S%qrdI, S%N)
            call permute_matrix_columns(S%qrdT, S%N, S%qrdP, S%qrdI)

            ! Q = full Q (from previous iteration QRP factorisation)
            call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

            do
                if (j .eq. S%L) then
                    exit
                else
                    ii = 0
                    do
                        if (ii .eq. S%north .or. j .eq. S%L) then
                            ! Q = Q * D (D from previous iteration)
                            call right_diagmult(S%qrdQ, S%qrdD, S%N)

                            ! QRP factorise Q
                            S%qrdP = 0
                            call dgeqp3(S%N, S%N, S%qrdQ, S%N, S%qrdP, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

                            ! D = diag(Q)
                            call diag(S%qrdQ, S%qrdD, S%N)

                            ! R = uppertri(Q)
                            call uppertri(S%qrdQ, S%qrdR, S%N)

                            ! R = inv(D) * R
                            call left_diaginvmult(S%qrdR, S%qrdD, S%N)

                            ! T = T * inv(P)
                            S%qrdI = S%qrdP
                            call invert_permutation(S%qrdI, S%N)
                            call permute_matrix_rows(S%qrdT, S%N, S%qrdI, S%qrdP)
            
                            ! T = R * T
                            call dtrmm('l', 'u', 'n', 'n', S%N, S%N, 1.0_dp, S%qrdR, S%N, S%qrdT, S%N)

                            ! Q = full Q (from previous iteration QRP factorisation)
                            call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

                            exit
                        else
                            j = j + 1
                            ii = ii + 1
                            call left_Bmult(S, S%qrdQ, getj(j, S%L, S%L), sigma)
                        endif
                    enddo
                endif
            enddo

            ! T = D * T
            call left_diagmult(S%qrdT, S%qrdD, S%N)
            ! R = Q * T
            call dgemm('n', 'n', S%N, S%N, S%N, 1.0_dp, S%qrdQ, S%N, S%qrdT, S%N, 0.0_dp, S%qrdR, S%N)   
            ! R = id + R
            do j = 1, S%N
                S%qrdR(j, j) = S%qrdR(j, j) + 1.0_dp
            enddo

            call determinant(Rold, S%qrdR, S%N, S%qrdP, info)


            S%h(i, l) = -S%h(i, l)



            ! Iteration j = 1
            j = 1

            ! Q = B(L)
            call make_B(S, S%qrdQ, getj(j, S%L, S%L), sigma)
                    
            if (j .lt. S%north) then
                do j = j+1, S%north
                    call left_Bmult(S, S%qrdQ, getj(j, S%L, S%L), sigma)
                enddo
            endif
            j = S%north

            ! QRP factorise Q
            S%qrdP = 0
            call dgeqp3(S%N, S%N, S%qrdQ, S%N, S%qrdP, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

            ! D = diag(Q)
            call diag(S%qrdQ, S%qrdD, S%N)

            ! T = uppertri(Q)
            call uppertri(S%qrdQ, S%qrdT, S%N)
            ! T = inv(D) * T
            call left_diaginvmult(S%qrdT, S%qrdD, S%N)
            ! T = T * inv(P)
            S%qrdI = S%qrdP
            call invert_permutation(S%qrdI, S%N)
            call permute_matrix_columns(S%qrdT, S%N, S%qrdP, S%qrdI)

            ! Q = full Q (from previous iteration QRP factorisation)
            call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

            do
                if (j .eq. S%L) then
                    exit
                else
                    ii = 0
                    do
                        if (ii .eq. S%north .or. j .eq. S%L) then
                            ! Q = Q * D (D from previous iteration)
                            call right_diagmult(S%qrdQ, S%qrdD, S%N)

                            ! QRP factorise Q
                            S%qrdP = 0
                            call dgeqp3(S%N, S%N, S%qrdQ, S%N, S%qrdP, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

                            ! D = diag(Q)
                            call diag(S%qrdQ, S%qrdD, S%N)

                            ! R = uppertri(Q)
                            call uppertri(S%qrdQ, S%qrdR, S%N)

                            ! R = inv(D) * R
                            call left_diaginvmult(S%qrdR, S%qrdD, S%N)

                            ! T = T * inv(P)
                            S%qrdI = S%qrdP
                            call invert_permutation(S%qrdI, S%N)
                            call permute_matrix_rows(S%qrdT, S%N, S%qrdI, S%qrdP)
            
                            ! T = R * T
                            call dtrmm('l', 'u', 'n', 'n', S%N, S%N, 1.0_dp, S%qrdR, S%N, S%qrdT, S%N)

                            ! Q = full Q (from previous iteration QRP factorisation)
                            call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

                            exit
                        else
                            j = j + 1
                            ii = ii + 1
                            call left_Bmult(S, S%qrdQ, getj(j, S%L, S%L), sigma)
                        endif
                    enddo
                endif
            enddo

            ! T = D * T
            call left_diagmult(S%qrdT, S%qrdD, S%N)
            ! R = Q * T
            call dgemm('n', 'n', S%N, S%N, S%N, 1.0_dp, S%qrdQ, S%N, S%qrdT, S%N, 0.0_dp, S%qrdR, S%N)   
            ! R = id + R
            do j = 1, S%N
                S%qrdR(j, j) = S%qrdR(j, j) + 1.0_dp
            enddo

            call determinant(Rnew, S%qrdR, S%N, S%qrdP, info)


            S%h(i, l) = -S%h(i, l)


            if (sigma .eq. 1) then      ! up spin (sigma = 1)
                S%deltaup = exp(-2 * sigma * S%alpha * S%h(i, l)) - 1
                S%Rup     = Rnew / Rold
            else                        ! dn spin (sigma = -1)
                S%deltadn = exp(-2 * sigma * S%alpha * S%h(i, l)) - 1
                S%Rdn     = Rnew / Rold
            endif


        endsubroutine scratch_R


endmodule metropolisratios_mod