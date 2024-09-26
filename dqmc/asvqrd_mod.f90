module asvqrd_mod
    use numbertypes
    use customla_mod
    implicit none

    ! TODO:
    ! Figure out how to pass subroutines as arguments, so the asvqrd
    ! algorithm can just be called with workspaces and a subroutine that
    ! returns the jth matrix in the chain (ie, a getB subroutine, see the
    ! implementation suggestion in subroutine chainiteration)

    contains


        subroutine chainiteration(B, Q, tau, R, D, T, I, P, work, lwork, info, N)
            !
            ! See the paper:
            !
            ! Stable solutions of linear systems involving long chain of matrix multiplications
            ! by: Zhaojun Bai, Che-Run Lee, Ren-Cang Li, Shufang Xu
            ! doi:10.1016/j.laa.2010.06.023
            !
            ! for a precise description of the chain multiplication algorithm.
            ! This subroutine implements lines 4 - 6 of the ASvQRD algorithm,
            ! which is a single iteration of the main loop of that algorithm. 
            !
            ! In the algorithm, those lines read (on iteration j, j-1 being the previous iteration):
            !
            ! C(j) = (B(j) * Q(j-1)) * D(j-1)     respect parenthesis
            ! C(j) = Q(j) * R(j) * P(j)           QR factorisation with column pivoting
            ! D(j) = diag(R(j))
            ! T(j) = inv(D(j)) * R(j) * P(j)
            !
            ! To save on memory and because of BLAS/LAPACK stucture, this iteration is
            ! carried out slightly differently than as stated in the paper.
            !
            ! Since (equation 2.5) the multiplication factors:
            !
            ! B(1) * ... * B(L) = Q(L) * D(L) * (T(L) * ... * T(1))
            !
            ! and the product of T matrices must be computed eventually (line 9 of ASvQRD),
            ! this subroutine updates:
            !
            ! T(j) = inv(D(j)) * R(j) * P(j) * T(j-1)
            !
            ! In other words, line 9 of ASvQRD is done on the fly.
            !
            ! In LAPACK, the QR factorisation with column pivoting (QRP factorisation)
            ! is done in the following form:
            !
            ! A * P = Q * R
            !
            ! or:
            !
            ! A = Q * R * inv(P)
            !
            ! So wherever there is multiplication by a permutation matrix P in ASvQRD, this subroutine will
            ! perform multiplication by the inverse of the calculated permutation matrix.
            !
            ! The LAPACK QRP subroutine dgeqp3 most notably acts as follows:
            !
            ! call dgeqp(A, P, tau)
            !
            ! In :
            !       A     : matrix to be QRP factorised
            !       P, tau: nothing needed on input
            !
            ! Out:
            !       A: overwritten, contains the following information:
            !          R         = uppertri(A)
            !          part of Q = strictlylowertri(A)
            !       P: contains information of A = Q * P * inv(P) permutation matrix
            !     tau: part of Q
            !
            ! Together, strictlylowertri(A) and tau contain the information of Q (eg, for explicitly
            ! constructing Q (dorgqr) or multiplying by Q (dormqr)).
            !
            ! How BLAS/LAPACK overwrites matrices turns out to be a headache for implementing/using
            ! the ASvQRD algorithm efficently. The following description is how I have implemented
            ! the algorithm (at least for the time being, better memory management is likely possible):
            !
            ! call dormqr(B, Q, tau)     (B = B * Q(j-1))
            ! B = B * D                  (B = B * D(j-1))
            ! call dgeqp3(B, P, tau)     (QRP factorise B)
            ! D = diag(B)                (D(j) = diag(R(j)))
            ! call uppertri(R, A)        (R = uppertri(A); explicit R with zeros below diagonal is needed for permutation multiplication)
            ! R = R * inv(P)             (T(j) = inv(D(j)) * R(j) * P(j) * T(j-1) steps)
            ! R = inv(D) * R
            ! T = R * T
            !
            ! Note that the Q(j) information is now stored in B and tau, and are needed
            ! as the Q and tau inputs in the j+1th iteration.
            ! This causes an implementation headache, where to iterate over this subroutine
            ! two matrices B and Q swap around:
            !
            ! iteration j  : B = B input, Q = Q input
            ! iteration j+1: B = Q input, Q = B input
            !
            ! Note that before calling this subroutine, the 1st iteration must be done separately (ie, 
            ! this subroutine handles j = 2, ..., L but not j = 1)
            !
            ! A full implementation of ASvQRD might then look as follows:
            !
            ! B = getB(1)
            ! call dgeqp3(B, P, tau)     (B is now Q input for iteration j = 2)
            ! D = diag(B)
            ! call uppertri(R, B)        (R = uppertri(B) with zeros below diagonal)
            ! R = R * inv(P)
            ! T = inv(D) * R             (iteration j = 1 complete)
            !
            ! do j = 2, L
            !      if     (j is even) then
            !           Q = getB(j)
            !           call chainiteration(B = Q, Q = B, tau, D, T)
            !      elseif (j is odd ) then
            !           B = getB(j)
            !           call chainiteration(B = B, Q = Q, tau, D, T)
            !      endif
            ! enddo
            !
            ! Now the final calculation of:
            !
            ! B(L) * .... * B(1) = Q(L) * D(L) * (T(L) * ... * T(1))
            !
            ! Using the last computed Q(L), D(L), and storing of T(L) * ... * T(1) in T
            !
            ! T = D * T
            !
            ! if     (L is even) then
            !      call dormqr(T, Q, tau)
            ! elseif (L is odd ) then
            !      call dormqr(T, B, tau)
            ! endif
            !
            ! Result: T = B(L) * ... * B(1)
            !
            ! Note:
            ! This implementation just considers forming the chain multiplication B(L) * ... * B(1) stabley,
            ! and not inversion.
            !
            ! A Green's function calculation requires the calculation of an inverse of the form:
            !
            ! G = inv(id + B(L) * ... * B(1))
            !
            ! This is calculated according to (see equation 2.11 of the paper):
            !
            ! G = inv((inv(Db) * trans(Q(L)) + Ds * T)) * inv(Db) * trans(Q(L))
            !
            ! where Db and Ds are diagonal matrices (stored as vectors) with values:
            !
            !         | D(i) if abs(D(i)) > 1
            ! Db(i) = | 1    otherwise
            !
            !         | D(i) if abs(D(i)) <= 1
            ! Ds(i) = | 1    otherwise
            !
            ! where D is the Lth calculated D matrix (ie, the last one).
            ! These matrices satisfy:
            !
            ! D = Db * Ds
            !
            ! A full implementation for computing a Green's function:
            !
            ! G = inv(id + B(L) * ... * B(1))
            !
            ! might then look as follows:
            !
            ! First calculating T such that Q(L) * D(L) * T = B(L) * ... * B(1)
            !
            ! B = getB(1)
            ! call dgeqp3(B, P, tau)     (B is now Q input for iteration j = 2)
            ! D = diag(B)
            ! call uppertri(R, B)        (R = uppertri(B) with zeros below diagonal)
            ! R = R * inv(P)
            ! T = inv(D) * R             (iteration j = 1 complete)
            !
            ! do j = 2, L
            !      if     (j is even) then
            !           Q = getB(j)
            !           call chainiteration(B = Q, Q = B, tau, D, T)
            !      elseif (j is odd ) then
            !           B = getB(j)
            !           call chainiteration(B = B, Q = Q, tau, D, T)
            !      endif
            ! enddo
            !
            ! Now computing the D = Db * Ds decomposition
            ! Storing Db in D and Ds in F
            !
            ! do j = 1, N
            !     if (abs(D(j)) .gt. 1) then
            !         F(j) = 1
            !     else
            !         F(j) = D(j)
            !         D(j) = 1
            !    endif
            ! enddo
            !
            ! Now to invert
            ! Trick: inv(Db) * trans(Q(L)) = trans(Q(L) * inv(Db))
            !
            ! if     (L is even) then            Q information stored in Q
            !            call dorgr(Q, tau)      Q = full Q
            !            Q = Q * inv(Db)         no transpose
            !            T = Ds * T              diagonal left multiply
            !            T = T + trans(Q)        double do loop
            !            T = inv(T)              PA = LU inversion
            !            G = T * trans(Q)        dgemm
            ! elseif (L is odd ) then            Q information stored in B
            !            call dorgr(B, tau)      B = full Q
            !            B = Q * inv(Db)         no transpose
            !            T = Ds * T              diagonal left multiply
            !            T = T + trans(Q)        double do loop
            !            T = inv(T)              PA = LU inversion
            !            G = T * trans(Q)        dgemm
            ! endif
            !
            real(dp), intent(inout) :: B(N, N)
            real(dp), intent(inout) :: Q(N, N)
            real(dp), intent(inout) :: tau(N)
            real(dp), intent(inout) :: R(N, N)
            real(dp), intent(inout) :: D(N)
            real(dp), intent(inout) :: T(N, N)
            integer , intent(inout) :: I(N)
            integer , intent(inout) :: P(N)
            real(dp), intent(inout) :: work(lwork)
            integer , intent(in)    :: lwork
            integer , intent(inout) :: info
            integer , intent(in)    :: N
            !
            ! On iteration j
            ! B:
            !      in: matrix to multiply by in chain
            !     out: part of Q information for  j+1 iteration (Q argument)
            ! Q:
            !      in: part of Q information from j-1 iteration
            !     out: destroyed
            ! tau:
            !      in: rest of Q information from j-1 iteration
            !     out: rest of Q information for  j+1 iteration (tau argument)
            ! R:
            !      in: workspace (for R part of upcoming QR factorisation)
            !     out: destroyed
            ! D:
            !      in: diag(R) from j-1 iteration
            !     out: diag(R) for  j+1 iteration (D argument)
            ! T:
            !      in: j-1th updated T
            !     out: jth   updated T (T argument)
            ! I:
            !      in: permutation workspace
            !     out: destroyed
            ! P
            !      in: permutation workspace
            !     out: destroyed
            ! work:
            !      in: workspace (for LAPACK)
            !     out: workspace
            ! lwork:
            !      in: length of work
            ! info:
            !      in: LAPACK info variable, nonzero entry can potentially cause issues
            !     out: destroyed
            ! N:
            !      in: size of square chain array
            !

            ! TODO:
            ! Figure out if there is a better way to manage memory
            ! See if more BLAS/LAPACK routines could be called/if they are better for multiplication (eg, 
            ! maybe use dlarscl2 for the update A = inv(D) * A)

            ! call dormqr(B, Q, tau)     (B = B * Q(j-1))
            call dormqr('R', 'N', N, N, N, Q, N, tau, B, N, work, lwork, info)

            ! B = B * D                  (B = B * D(j-1))
            call right_diagmult(B, D, N)

            ! call dgeqp3(B, P, tau)     (QRP factorise B)
            P = 0 ! zero pivots for LAPACK input
            call dgeqp3(N, N, B, N, P, tau, work, lwork, info)
            
            ! D = diag(B)                (D(j) = diag(R(j)))
            call diag(B, D, N)

            ! call uppertri(R, B)        (R = uppertri(B); explicit R with zeros below diagonal is needed for permutation multiplication)
            call uppertri(B, R, n)

            ! R = R * inv(P)             (T(j) = inv(D(j)) * R(j) * P(j) * T(j-1) steps)
            call invert_permutation(P, I, N) ! I = inv(P)
            call permutecols(R, I, N)        ! R = R * I

            ! R = inv(D) * R
            call left_diaginvmult(R, D, N)

            ! T = R * T
            call dlacpy('a', N, N, T, N, Q, N) ! Q = T, no A = B * A general matrix update subroutine in BLAS, copying T to Q to then use dgemm
            call dgemm('n', 'n', N, N, N, 1.0_dp, R, N, Q, N, 0.0_dp, T, N) ! T = R * Q
            ! It would be nice to find a way around this extra copying


        endsubroutine chainiteration


        subroutine firstchainiteration(B, tau, R, D, T, I, P, work, lwork, info, N)
            !
            ! Performs the first (j=1) iteration of the ASvQRD algorithm.
            ! See the subroutine chainiteration for more details.
            !
            ! Notably, the matrix B is overwritten, and becomes the Q input for the
            ! j=2 iteration.
            !
            real(dp), intent(inout) :: B(N, N)
            real(dp), intent(inout) :: tau(N)
            real(dp), intent(inout) :: R(N, N)
            real(dp), intent(inout) :: D(N)
            real(dp), intent(inout) :: T(N, N)
            integer , intent(inout) :: I(N)
            integer , intent(inout) :: P(N)
            real(dp), intent(inout) :: work(lwork)
            integer , intent(in)    :: lwork
            integer , intent(inout) :: info
            integer , intent(in)    :: N


            ! call dgeqp3(B, P, tau)     (B is now Q input for iteration j = 2)
            P = 0 ! zero pivots for LAPACK input
            call dgeqp3(N, N, B, N, P, tau, work, lwork, info)

            ! D = diag(B)
            call diag(B, D, N)

            ! call uppertri(R, B)        (R = uppertri(B) with zeros below diagonal)
            call uppertri(B, R, n)

            ! R = R * inv(P)
            call invert_permutation(P, I, N) ! I = inv(P)
            call permutecols(R, I, N)        ! R = R * I

            ! T = inv(D) * R             (iteration j = 1 complete)
            call left_diaginvmult(T, D, N)


        endsubroutine firstchainiteration


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


endmodule asvqrd_mod