    !>\brief Contains procedures for generating or updating the equal time greens functions \f$G_\sigma\f$.
    !!
    !! Only one of the following modules:
    !!
    !!     bmult_mod
    !!     bmultexact_mod
    !!
    !! should be used (ie, not commented out) in this module, depending on how multiplication
    !! by B matrices is desired to be done.
    !!
    !! The only difference in the two modules for multiplying by the \f$B_\sigma\f$ matrices is that
    !! bmult_mod uses the checkerboard method to approximate \f$\exp(\Delta\tau T)\f$ and \f$\exp(\Delta\tau T)^{-1}\f$,
    !! while bmultexact_mod uses no approximation (but runs much slower).
module equalgreens_mod
    use stduse
    use blas_interface
    use lapack_interface
    use simulationsetup_mod
    use bmult_mod
    use printing_mod
    use customla_mod
    use utilities
    ! use bmultexact_mod
    implicit none

    contains

        !> \brief Compares the currently held Green's function \f$G_\sigma\f$ with a newly computed one.
        !!
        !! Computes the difference between the stored Green's function \f$G_\sigma\f$.
        !! and a newly generated one for a given time slice `l` and spin `sigma`.
        !! Optionally replaces the stored Green's function with the new one based on the `keep` flag.
        !!
        !! \param[inout] S         (`Simulation`)                 Simulation data type.
        !! \param[in]    l         (`integer`)                    Time slice to compare at.
        !! \param[out]   diff      (`real(dp)`)                   Computed matrix difference between current \f$G_\sigma\f$ and a newly computed one.
        !! \param[in]    sigma     (`integer`)                    \f$\sigma\f$ in \f$G_\sigma\f$. `1` for \f$\uparrow\f$ and `-1` for \f$\downarrow\f$.
        !! \param[in]    keep      (`logical, Optional`)          If `.true.`, the simulation \f$G_\sigma\f$ is overwritten by the newly computed one and not if `.false.`. Default value: `.true.`.
        !! \param[in]    difftype  (`character(len=*), Optional`) Type of matrix difference to compute. See \ref matdiff for more details. Default value: nothing (matdiff default: `difflim`).
        subroutine compareG(S, l, diff, sigma, keep, difftype)
            type(simulation), intent(inout)           :: S
            integer         , intent(in)              :: l
            real(dp)        , intent(out)             :: diff
            integer         , intent(in)              :: sigma
            logical         , intent(in)   , optional :: keep
            character(len=*), intent(in)   , optional :: difftype

            logical  :: actualkeep

            if (.not. present(keep)) then
                actualkeep = .false.
            else
                actualkeep = keep
            endif

            if (sigma .eq. 1) then
                call copy_matrix(S%Gup, S%qrdQ, S%N)
                call newG(S, l, sigma)
                diff = matdiff(S%Gup, S%qrdQ, S%N, S%N, S%qrdB, difftype)
                if (.not. actualkeep) then
                    call copy_matrix(S%qrdQ, S%Gup, S%N)
                endif
            else
                call copy_matrix(S%Gdn, S%qrdQ, S%N)
                call newG(S, l, sigma)
                diff = matdiff(S%Gdn, S%qrdQ, S%N, S%N, S%qrdB, difftype)
                if (.not. actualkeep) then
                    call copy_matrix(S%qrdQ, S%Gdn, S%N)
                endif
            endif
        endsubroutine compareG

        subroutine flipupdate_alternative(S, i, sigma)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i
            integer         , intent(in)    :: sigma

            integer :: j, k

            if (sigma .eq. 1) then
                S%qrdB = S%Gup
                do j = 1, S%N
                    do k = 1, S%N
                        S%Gup(j, k) = S%qrdB(j, k) - (1.0_dp/S%Rup) * S%qrdB(j, i) * S%deltaup * (del(i, k) - S%qrdB(i, k))
                    enddo
                enddo
            else
                S%qrdB = S%Gdn
                do j = 1, S%N
                    do k = 1, S%N
                        S%Gdn(j, k) = S%qrdB(j, k) - (1.0_dp/S%Rdn) * S%qrdB(j, i) * S%deltadn * (del(i, k) - S%qrdB(i, k))
                    enddo
                enddo
            endif           
        endsubroutine flipupdate_alternative

        !> Updates \f$G_\sigma\f$ quickly if the Hubbard-Stratonovich field was flipped at site \f$i\f$ (assumed flipped and intermediate quantities computed, so no time slice \f$l\f$ needed here).
        !!
        !! Updates:
        !! \f[G_\sigma = G_\sigma + x\cdot\text{outer}(G_\sigma(:, i) - e_i, G(i, :))\f]
        !! where:
        !! \f[\begin{aligned}
        !! x             &= \frac{\delta_\sigma}{1 + \delta_\sigma(1-G_\sigma(i, i))}          \\
        !!               &= \frac{\delta_\sigma}{R_\sigma}                                     \\
        !! \delta_\sigma &= \exp(\sigma\alpha(h_{\text{new}}(i, l) - h_{\text{old}}(i, l)) - 1 \\
        !!               &= \exp(-2\sigma\alpha h_{\text{old}}(i, l)) - 1                      \\
        !! R_\sigma      &= 1 + \delta_\sigma(1 - G_\sigma(i, i))                              \\
        !! e_i           &= i\text{th Cartesian unit vector}                                   \\
        !! i             &= \text{site flipped}                                                \\
        !! l             &= \text{current imaginary time slice}
        !! \end{aligned}\f]
        !! The Hubbard-Stratonovich field \f$h\f$ is assumed to be flipped already.
        !! Specifically, if \f$h_{\text{old}}(i, l)\f$ is the old value of the Hubbard-Stratonovich
        !! field at site \f$i\f$ and imaginary time slice \f$l\f$ and \f$h_{\text{new}}(i, l)\f$ the new one,
        !! then \f$h(i, l) = h_{\text{new}}(i, l) = -h_{\text{old}}(i, l)\f$ (where \f$h(i, l)\f$ is the value
        !! of the Hubbard-Stratonovich field stored in the program as it is being run before calling this subroutine).
        !!
        !! The quick update technique technically requires both the site \f$i\f$ and imaginary time slice \f$l\f$
        !! to be known to be carried out, but it is assumed that the Metropolis weight \f$R_\sigma\f$ and
        !! factor \f$\delta_\sigma\f$ are already computed before calling this subroutine (eliminating \f$l\f$ dependence,
        !! these are calculated in \ref greens_r called before this subroutine).
        !!
        !! \param[inout] S      (`Simulation`) Simulation data type.
        !! \param[in]    i      (`integer`)    site number that `h` was flipped at (already flipped).
        !! \param[in]    sigma  (`integer`)    \f$\sigma\f$ in \f$G_\sigma\f$. `1` for \f$\uparrow\f$ and `-1` for \f$\downarrow\f$.
        !! \see Todo: potentially implement delayed update or submatrix update.
        subroutine flipupdate(S, i, sigma)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: i
            integer         , intent(in)    :: sigma

            real(dp) :: x

            if (sigma .eq. 1) then      ! up spin (sigma = 1)
                x = S%deltaup / S%Rup

                ! flipwork(:, 1) = G(:, i)
                call dcopy(S%N, S%Gup(1, i), 1, S%flipwork(1, 1), 1)

                ! flipwork(:, 2) = G(i, :) - ei
                call dcopy(S%N, S%Gup(i, 1), S%N, S%flipwork(1, 2), 1)
                S%flipwork(i, 2) = S%flipwork(i, 2) - 1.0_dp

                ! Outer product
                call dger(S%N, S%N, x, S%flipwork(1, 1), 1, S%flipwork(1, 2), 1, S%Gup, S%N)
            else                        ! dn spin (sigma = -1)
                x = S%deltadn / S%Rdn

                ! flipwork(:, 1) = G(:, i)
                call dcopy(S%N, S%Gdn(1, i), 1, S%flipwork(1, 1), 1)

                ! flipwork(:, 2) = G(i, :) - ei
                call dcopy(S%N, S%Gdn(i, 1), S%N, S%flipwork(1, 2), 1)
                S%flipwork(i, 2) = S%flipwork(i, 2) - 1.0_dp

                ! Outer product
                call dger(S%N, S%N, x, S%flipwork(1, 1), 1, S%flipwork(1, 2), 1, S%Gdn, S%N)
            endif
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
                call right_Binvmult(S, S%Gup, l, sigma)
                call left_Bmult(S, S%Gup, l, sigma)
            else                        ! dn spin (sigma = -1)
                call right_Binvmult(S, S%Gdn, l, sigma)
                call left_Bmult(S, S%Gdn, l, sigma)
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

            ! real(dp) :: diff
            ! real(dp) :: wrapG(S%N, S%N)

            ! Two choices for updating Green's function
            ! 1. From scratch (slow)
            ! call newG(S, l, sigma)
            ! return

            ! 2. Wrapping (fast)
            if (sigma .eq. 1) then          ! up spin (sigma = 1)
                if ((mod(S%upstabi+1, S%nstab) .eq. 0) .or. (l .eq. 1)) then
                    ! Todo: finish implementing wrap and scratch test
                    ! 
                    ! First wrap
                    ! call wrap(S, l, sigma)
                    ! wrapG = S%Gup

                    ! Now make G from scratch and see how well the wrapping did
                    ! call compareG(S, l, diff, sigma, keep=.true.)

                    ! Present for debugging ---------------------------------------
                    ! write(stdout, "(a, f17.8)") "Gup newG and wrap diff = ", diff
                    ! call print_matrix(S%Gup - wrapG, stdout, "Gup - wrap G = ")
                    ! call print_matrix(S%Gup, stdout, "Gup = ")
                    ! call print_matrix(wrapG, stdout, "wrap Gup = ")
                    ! -------------------------------------------------------------
                    
                    ! Only make a G from scratch
                    ! (no checking how wrapping did)
                    call newG(S, l, sigma)

                    S%upstabi = 0
                else
                    call wrap(S, l, sigma)
                    S%upstabi = S%upstabi + 1
                endif
            else                            ! dn spin (sigma = -1)
                if ((mod(S%dnstabi+1, S%nstab) .eq. 0) .or. (l .eq. 1)) then
                    ! Todo: finish implementing wrap and scratch test
                    ! 
                    ! First wrap
                    ! call wrap(S, l, sigma)
                    ! wrapG = S%Gdn

                    ! Now make G from scratch and see how well the wrapping did
                    ! call compareG(S, l, diff, sigma, keep=.true.)

                    ! Present for debugging ---------------------------------------
                    ! write(stdout, "(a, f17.8)") "Gdn newG and wrap diff = ", diff
                    ! call print_matrix(S%Gdn - wrapG, stdout, "Gdn - wrap G = ")
                    ! call print_matrix(S%Gdn, stdout, "Gdn = ")
                    ! call print_matrix(wrapG, stdout, "wrap Gdn = ")
                    ! -------------------------------------------------------------
                    
                    ! Only make a G from scratch
                    ! (no checking how wrapping did)
                    call newG(S, l, sigma)

                    S%dnstabi = 0
                else
                    call wrap(S, l, sigma)
                    S%dnstabi = S%dnstabi + 1
                endif
            endif
        endsubroutine timeupdate

        !> Computes a new equal time Green's function \f$G_\sigma(l)\f$ at time step \f$l\f$.
        !!
        !! Mathematically, the equal time Green's function \f$G_\sigma(l)\f$ at imaginary time step \f$l\f$ is:
        !! \f[G_\sigma(l) = (\text{id} + B_\sigma(l)\dots B_\sigma(1)B_\sigma(L)\dots B_\sigma(l+1))^{-1}\f]
        !! Naively done (ie, straight multiplication, addition of identity, inversion), this computation is unstable.
        !! Care must be taken to compute \f$G_\sigma(l)\f$ in a stable manner.
        !! This is done by the ASvQRD algorithm.
        !!
        !! The ASvQRD algorithm is as follows.
        !! To compute:
        !! \f[G = (\text{id} + B_L\dots B_1)^{-1}\f]
        !! (now written more generically), first the \f$B_j\f$ matrices are multiplied together stably
        !! by QRP factorizing after each multiplication from right to left, then the inverse is computed
        !! by taking advantage of a factored form of the \f$B_j\f$ matrices.
        !!
        !! Here is the algorithm (written in a way slightly modified from the original paper
        !! in a way better suited for code; the formatting is bad because I don't know how
        !! to format it well in Doxygen).
        !!
        !! \f$B_1P = QR\f$ </p>
        !! \f$D = \text{diag}(R)\f$ </p>
        !! \f$T = D^{-1}RP\quad(QRP)\f$ </p>
        !! \f$\text{do }j = 2,3,\dots,L:\f$ </p>
        !! \f$C = (B_jQ)D\f$ </p>
        !! \f$CP = QR\quad(QRP)\f$ </p>
        !! \f$D = \text{diag}(R)\f$ </p>
        !! \f$T = D^{-1}RP^{-1}T\f$ </p>
        !! \f$\text{enddo}\f$ </p>
        !! \f$D = D_bD_s\quad(\text{decompose})\f$ </p>
        !! \f$G = (D_sT + D_b^{-1}Q^T)^{-1}D_b^{-1}Q^T\f$
        !!
        !! The decomposition \f$D=D_bD_s\f$ is explained in the subroutine \ref dbds.
        !! After the main do loop iteration, the product of \f$B_j\f$ matrices is \f$B_L\dots B_1 = QDT\f$.
        !! A decomposition of \f$D\f$ and trick for inverting the sum is then used.
        !!
        !! The primary bottleneck of this subroutine (and DQMC itself) is the constant need
        !! to \f$QRP\f$ factorize. This is slightly alleviated through use of the `S%north`
        !! attribute of the `Simulation` data structure. `S%north` \f$B_j\f$ matrices are
        !! multiplied together before doing a \f$QRP\f$ factorisation (or one is done at
        !! the end if not enough have been multiplied). Care should be taken to not
        !! set this parameter too low and incur slowdown, or too high and incur instability.
        !!
        !! The \ref getj subroutine is used to ensure that the \f$B_\sigma(j)\f$ matrices
        !! are multiplied in the correct order to compute the Green's function \f$G_\sigma(l)\f$
        !! correctly.
        subroutine newG(S, l, sigma)
            type(Simulation), intent(inout) :: S
            integer         , intent(in)    :: l
            integer         , intent(in)    :: sigma

            integer :: j ! B matrix counter
            integer :: i ! north counter
            integer :: info

            ! Reminder about LAPACK's dgeqp3 --------------------------------------------
            ! (QR factorisation with column pivoting: QRP factorisation)
            ! 
            ! call dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
            !
            ! m    : rows in A
            ! n    : columns in A
            ! A    : matrix to factorize AP = QR
            ! lda  : leading dimension of A
            ! jpvt : permutation matrix P stored as a vector
            ! tau  : holds reflector information of Q
            ! work : workspace
            ! lwork: length of workspace
            ! info : information about call
            !
            ! After calling, A is overwritten.
            ! The upper triangular part of A contains R (R is an upper triangular matrix)
            ! The strictly lower triangular part of A, along with tau, contains
            ! information to use or construct in full Q.
            ! The permutation matrix P is stored as a vector jpvt.
            ! ---------------------------------------------------------------------------

            ! First iteration ---------------------------------------
            ! (needs some initial setup)

            ! Iteration j = 1
            j = 1

            ! Construct first B matrix
            ! (stored in Q)
            call make_B(S, S%qrdQ, getj(j, S%L, l), sigma)

            ! Multiply a total of north B matrices together
            ! (stored in Q)
            if (j .lt. S%north) then
                do j = j+1, S%north
                    call left_Bmult(S, S%qrdQ, getj(j, S%L, l), sigma)
                enddo
            endif
            j = S%north

            ! QRP factorise the product of B matrices
            ! (QRP factorise Q)
            S%qrdP = 0
            call dgeqp3(S%N, S%N, S%qrdQ, S%N, S%qrdP, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

            ! D = diag(Q) (diag(R) of QRP)
            call diag(S%qrdQ, S%qrdD, S%N)

            ! T = uppertri(Q) (R of QRP)
            call uppertri(S%qrdQ, S%qrdT, S%N)

            ! T = inv(D) * T  (inv(D) * R)
            call left_diaginvmult(S%qrdT, S%qrdD, S%N)

            ! T = T * inv(P)  ((inv(D) * R) * inv(P))
            S%qrdI = S%qrdP
            call invert_permutation(S%qrdI, S%N)
            call permute_matrix_columns(S%qrdT, S%N, S%qrdP, S%qrdI)

            ! Q = explicit Q (Q of QRP; from iteration j)
            call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

            ! End first iteration ----------------------------------------------------------

            ! Loop strategy:
            ! 
            ! Each time a B matrix is multiplied in, increment j and i.
            ! j counts how many B matrices have been multiplied in.
            ! i counts how many B matrices have been multiplied without doing a QRP factorisation.
            ! When i is north or j is L, QRP factorize and do the main ASvQRD step (j = L is special,
            ! to ensure a QRP factorisation is done before exiting the main loop).
            ! i is reset to 0 when this QRP factorisation is done.
            ! The outermost do loop tracks when j reaches L, when this happens we know both:
            !       a QRP factorisation was done in the last inner loop iteration
            !       L B matrices have been multiplied together.
            ! At this point the outermost do loop is exited.
            do
                if (j .eq. S%L) then
                    ! All B matrices multiplied and last result QRP factorised --> good to leave
                    exit
                else
                    ! A QRP factorisation will have been done before getting here.
                    ! Reset the north counter i to 0.
                    i = 0
                    do
                        ! Have north B matrices been multiplied together,
                        ! or we reached L B matrix multiplications and have to
                        ! clean up to leave?
                        if (i .eq. S%north .or. j .eq. S%L) then
                            ! Q = Q * D (D from previous QRP factorisation)
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

                            ! T = inv(P) * T
                            S%qrdI = S%qrdP
                            call invert_permutation(S%qrdI, S%N)
                            call permute_matrix_rows(S%qrdT, S%N, S%qrdI, S%qrdP)

                            ! T = R * T
                            call dtrmm('l', 'u', 'n', 'n', S%N, S%N, 1.0_dp, S%qrdR, S%N, S%qrdT, S%N)

                            ! Q = full Q (QRP factorisation done in this segment, for next time)
                            call dorgqr(S%N, S%N, S%N, S%qrdQ, S%N, S%qrdtau, S%qrdwork, S%qrdlwork, S%info)

                            ! Go back to outermost do loop
                            exit
                        else
                            ! Good to multiply in a B matrix without QRP factorizing.
                            j = j + 1
                            i = i + 1
                            call left_Bmult(S, S%qrdQ, getj(j, S%L, l), sigma)
                        endif
                    enddo
                endif
            enddo

            ! Inversion step

            ! D = Db * Ds decomposition
            ! D = Db, F = Ds
            call DbDs(S%qrdD, S%qrdF, S%N)

            ! Computing G = inv(Ds*T + inv(Db)*trans(Q))*inv(Db)*trans(Q)

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

        !>\brief Returns the index of the \f$i\f$th \f$B_\sigma\f$ matrix from the right in the
        !! multiplication chain when computed a single particle equal time Green's function \f$G_\sigma\f$
        !! from scratch.
        !!
        !! The \f$l\f$th imaginary time step Green's function \f$G_\sigma(l)\f$ is (note: the usual `L`
        !! instead of this subroutine's argument `LL` is being used here since two `L`/`l` arguments are
        !! not allowed):
        !! \f[G_\sigma(l) = (\text{id} + B_\sigma(l)\dots B_\sigma(1)B_\sigma(L)\dots B_\sigma(l+1))^{-1}\f]
        !! Computing this requires computing the matrix product \f$B_\sigma(l)\dots B_\sigma(1)B_\sigma(L)\dots B_\sigma(l+1)\f$
        !! in a stable way. This is done by the ASvQRD algorithm, (implemented in \ref newg) which builds up this
        !! product multiplying from right to left. To do this, \ref newg iterates an integer `i = 1, 2, ..., L`
        !! (`i` stands for the `i`th \f$B_\sigma\f$ matrix from the right). This subroutine ensures that the correct
        !! \f$B_\sigma\f$ matrix is chosen for a given `i`.
        !! \param[in] i  (`integer`) number matrix on the right (ie, `i = 1` is rightmost, `i = 2` second rightmost, ...).
        !! \param[in] LL (`integer`) \f$L\f$ in \f$G_\sigma(l)\f$ definition (number of imaginary time slices).
        !! \param[in] l  (`integer`) \f$l\f$ in \f$G_\sigma(l)\f$ definition (imaginary time slice under consideration).
        !! \return    j  (`integer`) index of \f$B_\sigma\f$ matrix corresponding to `i` (eg, when `i = 1`, `j = l+1`).
        function getj(i, LL, l) result(j)
            integer, intent(in)  :: i
            integer, intent(in)  :: LL
            integer, intent(in)  :: l

            integer :: j

            if (i .le. LL - l) then
                j = l + i
            else
                j = l - LL + i
            endif
        endfunction getj

        !> \brief Decomposes the final diagonal matrix from the ASvQRD algorithm as required.
        !!
        !! The ASvQRD algorithm asks for the last diagonal \f$D\f$ matrix (stored as a vector)
        !! to be decomposed as follows:
        !! \f[D_b(i) = \begin{cases}D(i) &\text{ if }|D(i)| > 1 \\ 1 &\text{ otherwise}\end{cases}\f]
        !! \f[D_s(i) = \begin{cases}D(i) &\text{ if }|D(i)| \le 1 \\ 1 &\text{ otherwise}\end{cases}\f]
        !! \param[inout] D (`real(dp), dimension(n)`) In: last diagonal \f$D\f$ matrix of ASvQRD stored as a vector. Out: diagonal \f$D_b\f$ matrix stored as a vector.
        !! \param[out]   F (`real(dp), dimension(n)`) Out: diagonal \f$D_s\f$ matrix stored as a vector.
        !! \param[in]    n (`integer`) Dimension of diagonal matrices involved.
        !! \see The ASvQRD algorithm is presented in \link https://www.sciencedirect.com/science/article/pii/S0024379510003198 \endlink.
        !! Equations 2.9 and 2.10 are relevant for this subroutine.
        subroutine DbDs(D, F, n)
            real(dp), intent(inout) :: D(n)
            real(dp), intent(out)   :: F(n)
            integer , intent(in)    :: n

            integer :: i

            do i = 1, N
                if (abs(D(i)) .gt. 1) then
                    F(i) = 1.0_dp
                else
                    F(i) = D(i)
                    D(i) = 1.0_dp
                endif
            enddo
        endsubroutine DbDs
endmodule equalgreens_mod