module equalgreens_modfast
    use numbertypes
    use simulationsetup_mod
    use bmult_mod
    ! use bmultexact_mod
    implicit none
    !
    ! Contains procedures for generating or updating the equal time greens functions Gup and Gdn.
    !
    ! Only one of the following modules:
    !
    !     bmult_mod
    !     bmultexact_mod
    !
    ! should be used (ie, not commented out) in this module, depending on how multiplication
    ! by B matrices is desired to be done.
    !
    ! The only difference in the two modules for multiplying by B matrices is that
    ! bmult_mod uses the checkerboard method to approximate exp(dtau * T) and inv(dtau * T),
    ! while bmultexact_mod uses no approximation (but runs much slower).
    !
    contains


        subroutine flipupdate(S, i, sigma)
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
            ! The Hubbard Stratonovich field h is assumed to be flipped
            ! at position (i, l) already. That is:
            !
            ! hnew(i, l) =  h(i, l)
            ! hold(i, l) = -h(i, l)
            !
            ! The ratio of weights R(up/dn) and factor delta(up/dn)
            ! are assumed to be calculated beforehand and stored in S
            !
            ! This is the rank k update algorithm (fast update algorithm) 
            ! for flipping a spin in the Hubbard model. In this DQMC,
            ! flipping h(i, l) --> -h(i, l) only has 1 nonzero entry in
            ! the delta updating matrix, so k = 1.
            !
            ! Mathematically, G is updated according to:
            !
            ! G = G + U * S * V
            !
            ! where:
            !
            ! U = G(:, i)                                                                    N   x (k=1)
            ! S = delta(i, i) * inv(T)  (k=1) x (k=1)  * (k=1) x (k=1)                 --> (k=1) x (k=1)
            ! T = id + Y * delta(i, i)  (k=1) x (k=1)  + (k=1) x (k=1) * (k=1) x (k=1) --> (k=1) x (k=1)
            ! V = (G(:, i) - ei) ** T      (N x (k=1)) **  T                           --> (k=1) x   N
            !
            ! Expanded out:
            !
            ! U = G(:, i)
            ! S = delta(i, i) * [1 - delta(i, i) * (G(i, i) - 1)]^-1
            !   = delta(i, i) / [1 - delta(i, i) * (G(i, i) - 1)]
            ! V = (G(i, :) - ei) ** T
            !
            ! Multiplied together:
            !
            ! U * S * V = [delta(i, i) / [1 - delta(i, i) * (G(i, i) - 1)]] * G(:, i) * (G(i, :) - ei) **T
            !           = x * outer(G(:, i), G(i, :) - ei)
            !
            ! where:
            !
            ! x = delta(i, i) / [1 - delta(i, i) * (G(i, i) - 1)]
            !   = delta(i, i) / [1 + (1 - G(i, i)) * delta(i, i)]
            !   = delta       / R
            !
            ! where, for this Hubbard model:
            !
            ! R     = 1 + (1 - G(i, i)) * delta
            ! delta = delta(i, i) = exp(-2 * sigma * alpha * h(i, l)) - 1
            !
            ! R being the Metropolis weight of accepting the flip h(i, l) --> -h(i, l)
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i
            integer         , intent(in)    :: sigma

            real(dp) :: x

            ! TODO:
            ! Implement delayed update


            if (sigma .eq. 1) then      ! up spin (sigma = 1)
                x = S%deltaup / S%Rup

                ! flipwork(:, 1) = G(:, i) - 1
                call dcopy(S%N, S%Gup(1, i), 1, S%flipwork(1, 1), 1)
                S%flipwork(i, 1) = S%flipwork(i, 1) - 1.0_dp

                ! flipwork(:, 2) = G(i, :)
                call dcopy(S%N, S%Gup(i, 1), S%N, S%flipwork(1, 2), 1)

                ! Outer product
                call dger(S%N, S%N, x, S%flipwork(1, 1), 1, S%flipwork(1, 2), 1, S%Gup, S%N)
            else                        ! dn spin (sigma = -1)
                x = S%deltadn / S%Rdn

                ! flipwork(:, 1) = G(:, i) - 1
                call dcopy(S%N, S%Gdn(1, i), 1, S%flipwork(1, 1), 1)
                S%flipwork(i, 1) = S%flipwork(i, 1) - 1.0_dp

                ! flipwork(:, 2) = G(i, :)
                call dcopy(S%N, S%Gdn(i, 1), S%N, S%flipwork(1, 2), 1)

                ! Outer product
                call dger(S%N, S%N, x, S%flipwork(1, 1), 1, S%flipwork(1, 2), 1, S%Gdn, S%N)
            endif


            ! No BLAS:

            ! integer :: j
            ! integer :: k

            ! if (sigma .eq. 1) then      ! up spin (sigma = 1)
                ! x = S%deltaup / S%Rup
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
            ! else                        ! dn spin (sigma = -1)
                ! x = S%deltadn / S%Rdn
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
            ! endif


        endsubroutine flipupdate


        subroutine wrap(S, l, sigma)
            !
            ! Updates:
            !
            ! G = B(l) * G * inv(B(l))
            !
            ! Where G is the up (sigma = 1) or  dn (sigma = -1) equal time Green's function
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
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            if (sigma .eq. 1) then      ! up spin (sigma = 1)
                ! Fast
                call right_Binvmult(S, S%Gup, l, sigma)
                call left_Bmult(S, S%Gup, l, sigma)
                ! Slow
                ! call right_Binvmultexact(S, S%Gup, l, sigma)
                ! call left_Bmultexact(S, S%Gup, l, sigma)
            else                        ! dn spin (sigma = -1)
                ! Fast
                call right_Binvmult(S, S%Gdn, l, sigma)
                call left_Bmult(S, S%Gdn, l, sigma)
                ! Slow
                ! call right_Binvmultexact(S, S%Gdn, l, sigma)
                ! call left_Bmultexact(S, S%Gdn, l, sigma)
            endif


        endsubroutine wrap


        subroutine timeupdate(S, l, sigma)
            !
            ! Updates G for a new time slice iteration in a sweep
            ! Where G is the up spin (sigma = 1) or dn spin (sigma) equal time Green's function
            !
            ! Updates G in one of the two ways:
            !
            ! if this is the nstabth time Gup is being updated since
            ! a new computation from scratch, then Gup is computed from scratch
            !
            ! otherwise Gup is updated from the previous time iteration by wrapping
            ! (even if the last iteration was in a separate sweep call)
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            ! Two choiced for updating Green's function
            ! 1. From scratch (slow)
            ! call newG(S, l, sigma)
            ! return

            ! 2. Wrapping (fast)
            if (sigma .eq. 1) then          ! up spin (sigma = 1)
                if ((mod(S%upstabi+1, S%nstab) .eq. 0) .or. (l .eq. 1)) then
                    call newG(S, l, sigma)
                    S%upstabi = 0
                else
                    call wrap(S, l, sigma)
                    S%upstabi = S%upstabi + 1
                endif
            else                            ! dn spin (sigma = -1)
                if ((mod(S%dnstabi+1, S%nstab) .eq. 0) .or. (l .eq. 1)) then
                    call newG(S, l, sigma)
                    S%dnstabi = 0
                else
                    call wrap(S, l, sigma)
                    S%dnstabi = S%dnstabi + 1
                endif
            endif


        endsubroutine timeupdate



        subroutine newG(S, l, sigma)
            !
            ! Computes a new equal time Green's function at imaginary time step l
            ! and sets it to S%Gup (sigma=1) or S%Gdn (sigma=-1)
            !
            ! Mathematically (ignoring spin):
            !
            !   G(l) = inv(id + B(l) * ... * B(1) * B(L) * ... * B(l+1))
            !
            !
            ! Where each B(j) is a single electron propogator from imaginary time step j-1 to j
            !
            ! This Green's function is computed stably by the ASvQRD algorithm.
            ! This algorithm consists of two steps for computing G(l).
            !
            ! First, the product of L matrices:
            !
            !       B(L) * B(L-1) * ... * B(1)
            !
            ! is computed in a stable manner.
            !
            ! Then the addition of the identity and inverse are taken at once,
            ! in a stable manner.
            !
            ! The product of L matrices is done in the following way.
            !
            ! First the matrix B(1) is QR with column pivoting factorized.
            !
            !       B(1) = Q * R * inv(P)
            !
            ! (an inverse of P is being used here to match LAPACK usage)
            !
            ! Implementation note: in DQMC we need to take a product of L B matrices.
            ! Each 
            !
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            integer :: j ! B matrix counter
            integer :: i ! north counter
            integer :: info


            ! Iteration j = 1
            j = 1

            ! Q = B(l)
            ! Fast
            call make_B(S, S%qrdQ, getj(j, S%L, l), sigma)
            ! Slow
            ! call make_Bexact(S, S%qrdQ, getj(j, S%L, l), sigma)
                    
            if (j .lt. S%north) then
                do j = j+1, S%north
                    ! Fast
                    call left_Bmult(S, S%qrdQ, getj(j, S%L, l), sigma)
                    ! Slow
                    ! call left_Bmultexact(S, S%qrdQ, getj(j, S%L, l), sigma)
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
                    i = 0
                    do
                        if (i .eq. S%north .or. j .eq. S%L) then
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
                            i = i + 1
                            ! Fast
                            call left_Bmult(S, S%qrdQ, getj(j, S%L, l), sigma)
                            ! Slow
                            ! call left_Bmultexact(S, S%qrdQ, getj(j, S%L, l), sigma)
                        endif
                    enddo
                endif
            enddo

            ! Now to finish calculating G

            ! D(L) = Db * Ds decomposition, see ASvQRD algorithm
            ! D = Db, F = Ds
            call DbDs(S%qrdD, S%qrdF, S%N)

            ! B = trans(Q)
             call trans(S%qrdB, S%qrdQ, S%N)
            ! B = inv(Db) * B
            call left_diaginvmult(S%qrdB, S%qrdD, S%N)

            ! T = Ds * T                  
            call left_diagmult(S%qrdT, S%qrdF, S%N)       
            ! T = T + B                                                
            call add_matrix(S%qrdT, S%qrdB, S%N)
            ! T = inv(T)
            call invert(S%qrdT, S%N, S%invP, S%invwork, S%invlwork, S%info)  

            ! G = T * B
            if (sigma .eq. 1) then ! sigma =  1 --> Gup
                call dgemm('n', 'n', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdB, S%N, 0.0_dp, S%Gup, S%N)   
            else                   ! sigma = -1 --> Gdn
                call dgemm('n', 'n', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdB, S%N, 0.0_dp, S%Gdn, S%N)
            endif

        endsubroutine newG












        ! Up for deletion
        subroutine newG_old(S, l, sigma)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

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
        
            ! Now to finish calculating G

            ! D(L) = Db * Ds decomposition, see ASvQRD algorithm
            ! D = Db, F = Ds
            call DbDs(S%qrdD, S%qrdF, S%N)

            ! Q = full Q
            call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)
            ! B = trans(Q)
            call trans(S%qrdB, S%qrdQ, S%N)
            ! B = inv(Db) * B
            call left_diaginvmult(S%qrdB, S%qrdD, S%N)
            ! T = Ds * T                  
            call left_diagmult(S%qrdT, S%qrdF, S%N)    
            ! T = T + B                                                
            call add_matrix(S%qrdT, S%qrdB, S%N)
            ! T = inv(T)
            call invert(S%qrdT, S%N, S%invP, S%invwork, S%invlwork, S%info)                             
            ! G = T * B
            if (sigma .eq. 1) then ! sigma =  1 --> Gup
                call dgemm('n', 'n', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdB, S%N, 0.0_dp, S%Gup, S%N)   
            else                   ! sigma = -1 --> Gdn
                call dgemm('n', 'n', S%N, S%N, S%N, 1.0_dp, S%qrdT, S%N, S%qrdB, S%N, 0.0_dp, S%Gdn, S%N)
            endif


        endsubroutine newG_old


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


endmodule equalgreens_modfast