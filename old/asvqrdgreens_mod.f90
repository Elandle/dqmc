module asvqrdgreens_mod
    use numbertypes
    use customla_mod
    use simulationsetup_mod
    use equalgreens_mod
    implicit none


    contains

        subroutine newG_stable(S, l, sigma)
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
            integer         , intent(in)    :: sigma

            integer :: j

            ! Iteration j = 1
            j = 1
            ! S%qrdQ = B(l)
            call make_B(S, S%qrdQ, getj(j, S%L, l), sigma)
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
            call invert_permutation(S%qrdP, S%qrdI, S%N) ! S%qrdI = inv(S%qrdP)
            call permutecols(S%qrdT, S%qrdI, S%N)
            
            do j = 2, S%L
                ! S%qrdQ = full Q (from previous iteration QRP factorisation)
                call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)
                ! S%qrdQ = B(getj(j)) * S%qrdQ
                call left_Bmult(S, S%qrdQ, getj(j, S%L, l), sigma)
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
                call invert_permutation(S%qrdP, S%qrdI, S%N) ! S%qrdI = inv(S%qrdP)
                call permutecols(S%qrdR, S%qrdI, S%N)
                ! S%qrdT = S%qrdR * S%qrdT
                call left_matmul(S%qrdT, S%qrdR, S%N, S%qrdB) ! Should try to get rid of calling this subroutine
                                                              ! S%qrdB: a work array to hold a matrix copy
            enddo
            
            ! Now S%qrdT = B(l) * ... * B(1) * B(L) * ... * B(l+1)
            ! Now to finish calculating G

            ! G = inv(id + qrdT)
            ! Q information stored in Q

            ! D(L) = Db * Ds decomposition, see ASvQRD algorithm
            ! D = Db, F = Ds
            call DbDs(S%qrdD, S%qrdF, S%N)

            ! Q = full Q
            call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)
            ! Q = Q * inv(Db), no transpose
            call right_diaginvmult(S%qrdQ, S%qrdD, S%N)               
            ! T = Ds * T                                  
            call left_diagmult(S%qrdT, S%qrdF, S%N)    
            ! T = T + trans(Q)                                                 
            call add_trans(S%qrdT, S%qrdQ, S%N)                         
            ! T = inv(T)                                
            call invert(S%qrdT, S%N, S%invP, S%invwork, S%invlwork, S%info)                             
            ! G = T * trans(Q)
            if (sigma .eq. 1) then ! sigma =  1 --> Gup
                call dgemm('n', 't', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdQ, S%N, 0.0_dp, S%Gup, S%N)   
            else                   ! sigma = -1 --> Gdn
                call dgemm('n', 't', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdQ, S%N, 0.0_dp, S%Gdn, S%N)
            endif




        endsubroutine newG_stable

        subroutine DbDs(D, F, n)
            !
            ! The ASvQRD algorithm asks for the last diagonal D matrix (stored as a vector)
            ! to be decomposed as follows:
            !
            !         | D(i) if abs(D(i)) > 1
            ! Db(i) = | 1    otherwise
            !
            !         | D(i) if abs(D(i)) <= 1
            ! Ds(i) = | 1    otherwise
            !
            ! This subroutine takes the last diagonal D matrix and sets it to Db,
            ! and takes another (unset) diagonal matrix F (stored as a vector) and sets it to Ds
            !
            real(dp), intent(inout) :: D(n)
            real(dp), intent(out)   :: F(n)
            integer , intent(in)    :: n

            integer :: i

            do i = 1, N
                if (abs(D(i)) .gt. 1) then
                    F(i) = 1
                else
                    F(i) = D(i)
                    D(i) = 1
                endif
            enddo


        endsubroutine DbDs


endmodule asvqrdgreens_mod