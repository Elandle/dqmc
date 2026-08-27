module convenientla_mod
    use stduse
    use customla_mod
    implicit none

    ! Contains convenient, not high-performing, linear algebra procedures.
    ! Mainly for testing
    !
    ! For example, the update:
    !
    !     A = A * B
    !
    ! without having to supply a workspace

    contains

        subroutine right_matrixmultiply(A, B, n)
            !
            ! Updates:
            !
            ! A = A * B
            !
            ! where A and B are n x n matrices
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: B(n, n)
            integer , intent(in)    :: n

            real(dp), allocatable :: C(:, :)

            allocate(C(n, n))

            ! C = A
            call copy_matrix(A, C, n)
            ! A = C * B
            call dgemm('n', 'n', n, n, n, 1.0_dp, C, n, B, n, 0.0_dp, A, n)
        endsubroutine

        subroutine left_matrixmultiply(A, B, n)
            !
            ! Updates:
            !
            ! A = B * A
            !
            ! where A and B are n x n matrices
            !
            real(dp), intent(inout) :: A(n, n)
            real(dp), intent(in)    :: B(n, n)
            integer , intent(in)    :: n

            real(dp), allocatable :: C(:, :)

            allocate(C(n, n))

            ! C = A
            call copy_matrix(A, C, n)
            ! A = B * C
            call dgemm('n', 'n', n, n, n, 1.0_dp, B, n, C, n, 0.0_dp, A, n)
        endsubroutine

        real(dp) function twonorm(A) result(norm)
            !
            ! Returns the 2 norm of the matrix A
            !
            real(dp), intent(in)  :: A(:, :)
            real(dp), allocatable :: B(:, :)
            real(dp), allocatable :: S(:)
            real(dp), allocatable :: work(:)
            integer               :: m
            integer               :: n
            integer               :: lwork
            integer               :: info
            
            m = size(A, 1)
            n = size(A, 2)
            lwork = 5 * (3 * min(m, n) + max(m, n) + 5 * min(m, n))

            allocate(B(m, n))
            allocate(S(min(m, n)))
            allocate(work(lwork))

            call dlacpy('A', m, n, A, m, B, m)
            call dgesvd('N', 'N', m, n, B, m, S, work, lwork, work, lwork, work, lwork, info)

            norm = S(1) 
        endfunction twonorm

        subroutine dgpadm(ideg, m, t, H, ldh, wsp, lwsp, ipiv, iexph, ns, iflag)
            integer  :: ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
            real(dp) :: t, H(ldh,m), wsp(lwsp)
            !
            ! Adapted from the expokit subroutine dgpadm
            !
            ! -----Purpose----------------------------------------------------------|
            !
            !      Computes exp(t*H), the matrix exponential of a general matrix in
            !      full, using the irreducible rational Pade approximation to the 
            !      exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
            !      combined with scaling-and-squaring.
            !
            ! -----Arguments--------------------------------------------------------|
            !
            !      ideg      : (input) the degre of the diagonal Pade to be used.
            !                 a value of 6 is generally satisfactory.
            !
            !      m         : (input) order of H.
            !
            !      H(ldh,m)  : (input) argument matrix.
            !
            !      t         : (input) time-scale (can be < 0).
            !                  
            !      wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ ideg+1.
            !
            !      ipiv(m)   : (workspace)
            !
            ! >>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
            !                  i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
            !                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !                  NOTE: if the routine was called with wsp(iptr), 
            !                        then exp(tH) will start at wsp(iptr+iexph-1).
            !
            !      ns        : (output) number of scaling-squaring used.
            !
            !      iflag     : (output) exit flag.
            !                       0 - no problem
            !                      <0 - problem
            !
            ! ----------------------------------------------------------------------|
            !      Roger B. Sidje (rbs@maths.uq.edu.au)
            !      EXPOKIT: Software Package for Computing Matrix Exponentials.
            !      ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
            ! ----------------------------------------------------------------------|
            !
            integer  :: mm, i, j, k, ih2, ip, iq, iused, ifree, iodd, icoef, iput, iget
            real(dp) :: hnorm, scale, scale2, cp, cq

            ! intrinsic INT,ABS,DBLE,LOG,MAX

            ! ---  check restrictions on input parameters ...
            mm = m * m
            iflag = 0
            if (ldh .lt. m) iflag = -1
            if (lwsp .lt. 4*mm + ideg+1 ) iflag = -2
            if (iflag .ne. 0 ) stop "bad sizes (in input of dgpadm)"
            !
            ! ---  initialise pointers ...
            !
            icoef = 1
            ih2 = icoef + (ideg+1)
            ip  = ih2 + mm
            iq  = ip + mm
            ifree = iq + mm
            !
            ! ---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
            !      and set scale = t/2^ns ...
            !
            do i = 1, m
                wsp(i) = 0.0_dp
            enddo

            do j = 1, m
                do i = 1, m
                    wsp(i) = wsp(i) + abs(H(i, j))
                enddo
            enddo

            hnorm = 0.0_dp
            do i = 1, m
                hnorm = max(hnorm, wsp(i))
            enddo
            hnorm = abs(t * hnorm)

            if (hnorm .eq. 0.0_dp) stop "Error - null H in input of dgpadm."
            ns = max(0, int(log(hnorm)/log(2.0_dp)) + 2)
            scale = t / real(2 ** ns, dp)
            scale2 = scale * scale
            !
            ! ---  compute Pade coefficients ...
            !
            i = ideg + 1
            j = 2*ideg + 1
            wsp(icoef) = 1.0_dp
            do k = 1, ideg
                wsp(icoef+k) = (wsp(icoef+k-1) * real(i-k, dp)) / real(k*(j-k), dp)
            enddo
            !
            ! ---  H2 = scale2 * H * H ...
            !
            call dgemm('n', 'n', m, m, m, scale2, H, ldh, H, ldh, 0.0_dp, wsp(ih2), m)
            !
            ! ---  initialize p (numerator) and q (denominator) ...
            !
            cp = wsp(icoef+ideg-1)
            cq = wsp(icoef+ideg)
            do j = 1, m
                do i = 1, m
                    wsp(ip + (j-1)*m + i-1) = 0.0_dp
                    wsp(iq + (j-1)*m + i-1) = 0.0_dp
                enddo
                wsp(ip + (j-1)*(m+1)) = cp
                wsp(iq + (j-1)*(m+1)) = cq
            enddo
            !
            ! ---  Apply Horner rule ...
            !
            iodd = 1
            k = ideg - 1
100         continue
            iused = iodd*iq + (1-iodd)*ip
            call dgemm('n', 'n', m, m, m, 1.0_dp, wsp(iused), m, wsp(ih2), m, 0.0_dp, wsp(ifree), m)
            do j = 1, m
                wsp(ifree+(j-1)*(m+1)) = wsp(ifree + (j-1)*(m+1)) + wsp(icoef + k-1)
            enddo
            ip = (1-iodd)*ifree + iodd*ip
            iq = iodd*ifree + (1-iodd)*iq
            ifree = iused
            iodd = 1-iodd
            k = k-1
            if (k .gt. 0)  goto 100
            !
            ! ---  Obtain (+/-)(I + 2*(p\q)) ...
            !
            if (iodd .eq. 1) then
                call dgemm('n', 'n', m, m, m, scale, wsp(iq), m, H, ldh, 0.0_dp, wsp(ifree), m)
                iq = ifree
            else
                call dgemm('n', 'n', m, m, m, scale, wsp(ip), m, H, ldh, 0.0_dp, wsp(ifree), m)
                ip = ifree
            endif
            call daxpy(mm, -1.0_dp, wsp(ip), 1, wsp(iq), 1)
            call dgesv(m, m, wsp(iq), m, ipiv, wsp(ip), m, iflag)
            if (iflag .ne. 0) stop "Problem in dgesv (within dgpadm)"
            call dscal(mm, 2.0_dp, wsp(ip), 1)
            do j = 1, m
                wsp(ip + (j-1)*(m+1)) = wsp(ip + (j-1)*(m+1)) + 1.0_dp
            enddo
            iput = ip
            if (ns .eq. 0 .and. iodd .eq.1) then
                call dscal(mm, -1.0_dp, wsp(ip), 1)
                goto 200
            endif
            !
            ! --   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
            !
            iodd = 1
            do k = 1, ns
                iget = iodd*ip + (1-iodd)*iq
                iput = (1-iodd)*ip + iodd*iq
                call dgemm('n', 'n', m, m, m, 1.0_dp, wsp(iget), m, wsp(iget), m, 0.0_dp, wsp(iput), m)
                iodd = 1 - iodd
            enddo
200         continue
            iexph = iput
        endsubroutine dgpadm
        ! ----------------------------------------------------------------------|

        subroutine expm(A, exptA, t, ideg)
            !
            ! Convenient form of expokit's dgpadm (not good for heavy use)
            ! Sets exptA = exp(t*A) with A unmodified using a degree ideg
            ! Pade approximant in expokit's dgpadm.
            ! t and ideg are optional. By default, t = 1 and ideg = 6
            ! (a typical one suggested by expokit).
            !
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: exptA(size(A, 2), size(A, 2))
            real(dp), optional    :: t
            integer , optional    :: ideg

            integer  :: m
            integer  :: lda
            integer  :: lwsp
            integer  :: iflag
            integer  :: ns
            integer  :: iexpa
            integer  :: ipiv(size(A, 2))
            real(dp) :: ta
            integer  :: idega
            real(dp), allocatable :: wsp(:)

            if (present(t)) then
                ta = t
            else
                ta = 1.0_dp
            endif

            if (present(ideg)) then
                idega = ideg
            else
                idega = 6
            endif
            
            lda = size(A, 1)
            m = size(A, 2)
            lwsp = 4*m*m + idega+1    + 4*m*m
            iflag = 0
            ns = 0
            allocate(wsp(lwsp))

            call dgpadm(idega, m, ta, A, lda, wsp, lwsp, ipiv, iexpa, ns, iflag)
            call dcopy(m*m, wsp(iexpa), 1, exptA, 1)

            deallocate(wsp)
        endsubroutine expm

        function inverse(A) result(invA)
            real(dp), intent(in)  :: A(:, :)

            real(dp) :: invA(size(A, 1), size(A, 2))
            integer  :: P(size(A, 1))
            real(dp) :: work(64 * size(A, 1))
            integer  :: info

            call dcopy(size(A, 1)*size(A, 1), A, 1, invA, 1)
            call dgetrf(size(invA, 1), size(invA, 1), invA, size(invA, 1), P, info)
            call dgetri(size(invA, 1), invA, size(invA, 1), P, work, 4*size(invA, 1), info)
        endfunction inverse
endmodule convenientla_mod