!> \brief Contains procedures for exactly computing the matrix exponential.
module expm_mod
    use iso_fortran_env, only: real64
    implicit none

    integer, parameter, private :: dp = real64

    contains

        !> Adapted from the expokit subroutine dgpadm
        !!
        !! -----Purpose----------------------------------------------------------|
        !!
        !!      Computes exp(t*H), the matrix exponential of a general matrix in
        !!      full, using the irreducible rational Pade approximation to the 
        !!      exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
        !!      combined with scaling-and-squaring.
        !!
        !! -----Arguments--------------------------------------------------------|
        !!
        !!      ideg           : (input) the degre of the diagonal Pade to be used.
        !!                               a value of 6 is generally satisfactory.
        !!
        !!      m              : (input) order of H.
        !!
        !!      H(ldh,m)       : (input) argument matrix.
        !!
        !!      t              : (input) time-scale (can be < 0).
        !!                  
        !!      wsp(lwsp)      : (workspace/output) lwsp .ge. 4*m*m+ ideg+1.
        !!
        !!      ipiv(m)        : (workspace)
        !!
        !!      >>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
        !!                                i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
        !!                                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !!                                NOTE: if the routine was called with wsp(iptr), 
        !!                                      then exp(tH) will start at wsp(iptr+iexph-1).
        !!
        !!      ns             : (output) number of scaling-squaring used.
        !!
        !!      iflag          : (output) exit flag.
        !!                                 0 - no problem
        !!                                <0 - problem
        !!
        !! ----------------------------------------------------------------------|
        !!
        !! Roger B. Sidje (rbs@maths.uq.edu.au)
        !!
        !! EXPOKIT: Software Package for Computing Matrix Exponentials.
        !!
        !! ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
        !!
        !! ----------------------------------------------------------------------|
        subroutine dgpadm(ideg, m, t, H, ldh, wsp, lwsp, ipiv, iexph, ns, iflag)
            integer  :: ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
            real(dp) :: t, H(ldh,m), wsp(lwsp)
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


        !> Convenient form of EXPOKIT's `dgpadm` (not good for heavy use).
        !!
        !! Sets `exptA = exp(t*A)` with `A` unmodified using a degree `ideg`
        !! Pade approximant in EXPOKIT's `dgpadm`.
        !! `t` and `ideg` are optional.
        !! By default `t = 1` and `ideg = 6` (a typical `ideg` value suggest by EXPOKIT).
        !!
        !! \param[in]  A     (`real(dp), dimension(:, :)`)                   Matrix to diagonalize (should be square `size(A, 1) = size(A, 2)`).
        !! \param[out] exptA (`real(dp), dimension(size(A, 2), size(A, 2))`) Matrix to hold computed `exp(t*A)`.
        !! \param[in]  t     (`real(dp), optional`)                          Scalar `t` in `exp(t*A)`. Default value: `1`.
        !! \param[in]  ideg  (`integer, optional`)                           Degree of Pade approximant to use in `dgpadm`. Default value: `6`.
        subroutine expm(A, exptA, t, ideg)
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
endmodule expm_mod