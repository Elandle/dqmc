module polarization_mod
    use numbertypes
    use simulationsetup_mod
    use equalgreens_mod
    implicit none

    contains



        subroutine P_mult(A, B, P, N)
            !
            ! Sets:
            !
            !       A = B * P
            !
            ! where A is a complex matrix, B is a real matrix, and
            ! P is a complex diagonal matrix stored as a vector.
            !
            ! For use in computing B * P, where P is the matrix made by make_P,
            ! and B is the matrix:
            !
            !       B = B(L) * ... * B(1)
            !
            ! where each B(l) the lth single particle propogator of a particular spin
            !
            complex(dp), intent(out) :: A(N, N)
            real   (dp), intent(in)  :: B(N, N)
            complex(dp), intent(in)  :: P(N)
            integer    , intent(in)  :: N

            integer :: i

            ! A = G
            call zlacp2('a', N, N, B, N, A, N)

            ! A = A * P
            do i = 1, N
                call zscal(n, P(i), A(1, i), 1)
            enddo


        endsubroutine P_mult


        subroutine newB(S, l, sigma, B)
            !
            ! newB and newB2 do not agree.... needs fixing
            !
            !
            ! Sets:
            !
            !       B = B(l) * ... * B(1) * B(L) * ... * B(l+1)
            !
            ! where each B(l) is the lth single particle propogator of a particular
            ! spin sigma.
            !
            ! For use in constructing the matrix:
            !
            !       B = B(L) * ... * B(1)
            !
            ! needed in Polarization measurements.
            ! By determinant properties, the particular l input does not matter for
            ! these measurements, but to have all the B(l)'s in order, input:
            !
            !       l = L
            !
            ! where L is the number of imaginary time steps (stored in S%L).
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma
            real(dp)        , intent(out)   :: B(S%N, S%N)

            integer :: j ! B matrix counter
            integer :: i ! north counter

            ! Iteration j = 1
            j = 1
            ! S%qrdQ = B(l)
            call make_B(S, S%qrdQ, getj(j, S%L, l), sigma)
            ! north-1 more B multiplications before QRP factorisation
            do j = 2, S%north
                call left_Bmult(S, S%qrdQ, getj(j, S%L, l), sigma)
            enddo
            ! QRP factorise S%qrdQ
            S%qrdP = 0
            call dgeqp3(S%N, S%N, S%qrdQ, S%N, S%qrdP, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)
            ! S%qrdD = diag(S%qrdQ)
            call diag(S%qrdQ, S%qrdD, S%N)
            ! S%qrdT = uppertri(S%qrdQ)
            call uppertri(S%qrdQ, S%qrdT, S%N)
            ! S%qrdT = inv(S%qrdD) * S%qrdT
            call left_diaginvmult(S%qrdT, S%qrdD, S%N)
            ! S%qrdT = S%qrdT * inv(P)
            call colpivswap(S%qrdT, S%qrdP, S%N, S%qrdB)

            ! north B matrices have been multiplied and QRP factorised by now
            j = S%north + 1
            do while (j .lt. S%L)
                ! S%qrdQ = full Q (from previous iteration QRP factorisation)
                call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)
                ! Multiply north B matrices into Q (unless j reaches L along the way)
                i = 0
                do while ((i .lt. S%north) .and. (j .lt. S%L))
                    ! S%qrdQ = B(getj(j)) * S%qrdQ
                    call left_Bmult(S, S%qrdQ, getj(j, S%L, l), sigma)
                    i = i + 1
                    j = j + 1
                enddo
                ! S%qrdQ = S%qrdQ * S%qrdD (S%qrdD from previous iteration)
                call right_diagmult(S%qrdQ, S%qrdD, S%N)
                ! QRP factorise S%qrdQ
                S%qrdP = 0
                call dgeqp3(S%N, S%N, S%qrdQ, S%N, S%qrdP, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)
                ! S%qrdD = diag(S%qrdQ)
                call diag(S%qrdQ, S%qrdD, S%N)
                ! S%qrdR = uppertri(S%qrdQ)
                call uppertri(S%qrdQ, S%qrdR, S%N)
                ! S%qrdR = inv(S%qrdD) * S%qrdR
                call left_diaginvmult(S%qrdR, S%qrdD, S%N)
                ! S%qrdR = S%qrdR * inv(P)
                call colpivswap(S%qrdR, S%qrdP, S%N, S%qrdmatwork)
                ! S%qrdT = S%qrdR * S%qrdT
                call left_matmul(S%qrdT, S%qrdR, S%N, S%qrdB) ! Should try to get rid of calling this subroutine
                                                              ! S%qrdB: a work array to hold a matrix copy
            enddo
            !
            ! Now:
            !
            !       B = Q * D * T
            !
            ! Now to carry out these multiplications to compute B
            !
            ! Note, from the previous iteration:
            !
            !       Q information stored in S%qrdtau and
            !         the strictly lower triangular part of S%qrdQ
            !       D (diagonal matrix) stored in the vector S%qrdD
            !       T stored in S%qrdT
            !

            !
            ! TODO: rewrite subroutine just for making B (not adapted),
            !       notably make S%qrdT be B instead (avoid a copy)
            !
            !       think about reordering matrix products
            !
            !       B is needed to compute:
            !
            !           det(id + B)
            !           det(id + B * P)
            !   
            !       think if these can be done faster or more stably
            !       
            
            ! T = D * T
            ! would it be better to multiply by R on the right instead of left?
            call left_diagmult(S%qrdT, S%qrdD, S%N)

            ! T = Q * T
            call dormqr('L', 'N', S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdT, S%N, S%qrdwork, S%qrdlwork, S%info)

            ! B = T
            call copy_matrix(S%qrdT, B, S%N)


        endsubroutine newB


        subroutine newB2(S, l, sigma, B)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma
            real(dp)        , intent(out)   :: B(S%N, S%N)

            integer :: j

            call make_identity(B, S%N)

            do j = 1, S%L
                call left_Bmult(S, B, getj(j, S%L, l), sigma)
            enddo


        endsubroutine newB2


        subroutine measure_P(S, B, BP, p, sigma)
            !
            ! Measures Polarization:
            !
            !       det(id + B(L) * ... * B(1) * exp(2*pi*i/L * P)) / det(id + B(L) * ... * B(1))
            !
            ! from scratch (to be made more efficient later once it is confirmed it works).
            !
            ! This is measured by:
            !
            !       computing B = B(L) * ... * B(1)
            !       computing det(id + B)
            !       computing det(id + B * exp(2*pi*i/L) * P)
            !
            ! with no shortcuts.
            !
            ! The measurement is stored in p.
            !
            type(Simulation), intent(inout) :: S
            real(dp)        , intent(out)   :: B(S%N, S%N)
            complex(dp)     , intent(out)   :: BP(S%N, S%N)
            complex(dp)     , intent(out)   :: p
            integer         , intent(in)    :: sigma

            integer     :: i
            real(dp)    :: detB
            complex(dp) :: detBP

            ! TODO:
            ! find a more stable way to compute det(id + B)
            ! than just plain adding id to B and computing the determinant
            ! by LU

            ! B = B(L) * ... * B(1)
            call newB(S, S%L, sigma, B)
            ! BP = B * exp(2*pi*i/L) * P
            call P_mult(BP, B, S%polP, S%N)

            ! B  = B  + id
            ! BP = BP + id
            do i = 1, S%N
                B (i, i) = B (i, i) + 1.0_dp
                BP(i, i) = BP(i, i) + 1.0_dp
            enddo

            call  determinant(detB , B , S%N, S%qrdP, S%info)
            call zdeterminant(detBP, BP, S%N, S%qrdP, S%info)
            
            p = detBP / detB

        endsubroutine measure_P



endmodule polarization_mod