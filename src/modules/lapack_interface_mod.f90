module lapack_interface
    use stduse
    implicit none

    ! -----------------------------------------------------------------------------------------------------
    ! Linear solve, AX = B --------------------------------------------------------------------------------

    interface
        subroutine dgetrf(m   , n, a, lda, ipiv,   &
                          info)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            real(dp), intent(inout) :: a(lda, *)
            integer , intent(in)    :: lda
            integer , intent(out)   :: ipiv(*)
            integer , intent(out)   :: info
        endsubroutine dgetrf
    endinterface

    interface
        subroutine zgetrf(m   , n, a, lda, ipiv,    &
                          info)
            import                     :: dp
            integer    , intent(in)    :: m
            integer    , intent(in)    :: n
            complex(dp), intent(inout) :: a(lda, *)
            integer    , intent(in)    :: lda
            integer    , intent(out)   :: ipiv(*)
            integer    , intent(out)   :: info
        endsubroutine zgetrf
    endinterface

    interface
        subroutine dgetri(n    , a   , lda, ipiv, work,    &
                          lwork, info)
            import                  :: dp
            integer , intent(in)    :: n
            real(dp), intent(inout) :: a(lda, *)
            integer , intent(in)    :: lda
            integer , intent(in)    :: ipiv(*)
            real(dp), intent(out)   :: work(*)
            integer , intent(in)    :: lwork
            integer , intent(out)   :: info
        endsubroutine dgetri
    endinterface

    ! end Linear solve, AX = B ----------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Least squares ---------------------------------------------------------------------------------------

    ! end Least squares -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Orthogonal/unitary factors (QR, CS, etc.) -----------------------------------------------------------

    interface
        subroutine dgeqp3(m  , n   , a    , lda , jpvt,     &
                          tau, work, lwork, info)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            real(dp), intent(inout) :: A(lda, *)
            integer , intent(in)    :: lda
            integer , intent(inout) :: jpvt(*)
            real(dp), intent(out)   :: tau(*)
            real(dp), intent(out)   :: work(*)
            integer , intent(in)    :: lwork
            integer , intent(out)   :: info
        endsubroutine dgeqp3
    endinterface

    interface
        subroutine dorgqr(m  , n   , k    , a   , lda,      &
                          tau, work, lwork, info)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            integer , intent(in)    :: k
            real(dp), intent(inout) :: A(lda, *)
            integer , intent(in)    :: lda
            real(dp), intent(in)    :: tau(*)
            real(dp), intent(out)   :: work(*)
            integer , intent(in)    :: lwork
            integer , intent(out)   :: info
        endsubroutine dorgqr
    endinterface

    ! end Orthogonal/unitary factors (QR, CS, etc.) -------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Non-symmetric eigenvalues ---------------------------------------------------------------------------

    interface
        subroutine dgeev(jobvl, jobvr, n    , a   , lda,    &
                         wr   , wi   , vl   , ldvl, vr ,    &
                         ldvr , work , lwork, info)
            import                   :: dp
            character, intent(in)    :: jobvl
            character, intent(in)    :: jobvr
            integer  , intent(in)    :: n
            real(dp) , intent(inout) :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(out)   :: wr(*)
            real(dp) , intent(out)   :: wi(*)
            real(dp) , intent(out)   :: vl(ldvl, *)
            integer  , intent(in)    :: ldvl
            real(dp) , intent(out)   :: vr(ldvr, *)
            integer  , intent(in)    :: ldvr
            real(dp) , intent(out)   :: work(*)
            integer  , intent(in)    :: lwork
            integer  , intent(out)   :: info
        endsubroutine dgeev
    endinterface

    ! end Non-symmetric eigenvalues -----------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Hermitian/symmetric eigenvalues ---------------------------------------------------------------------

    ! end Hermitian/symmetric eigenvalues -----------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Singular Value Decomposition (SVD) ------------------------------------------------------------------

    ! end Singular Value Decomposition (SVD) --------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! BLAS-like -------------------------------------------------------------------------------------------

    interface
        subroutine dlaset(uplo, m  , n, alpha, beta,  &
                          a   , lda)
            import                 :: dp
            character, intent(in)  :: uplo
            integer  , intent(in)  :: m
            integer  , intent(in)  :: n
            real(dp) , intent(in)  :: alpha
            real(dp) , intent(in)  :: beta
            real(dp) , intent(out) :: a(lda, *)
            integer  , intent(in)  :: lda
        endsubroutine dlaset
    endinterface

    interface
        subroutine dlascl2(m, n, d, x, ldx)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            integer , intent(in)    :: d(*)
            real(dp), intent(inout) :: x(ldx, *)
            integer , intent(in)    :: ldx
        endsubroutine dlascl2
    endinterface

    interface
        subroutine dlacpy(uplo, m  , n, a, lda,   &
                          b   , ldb)
            import                 :: dp
            character, intent(in)  :: uplo
            integer  , intent(in)  :: m
            integer  , intent(in)  :: n
            real(dp) , intent(in)  :: a(lda, *)
            integer  , intent(in)  :: lda
            real(dp) , intent(out) :: b(ldb, *)
            integer  , intent(in)  :: ldb
        endsubroutine dlacpy
    endinterface

    interface
        real(dp) function dlange(norm, m, n, a, lda,    &
                                 work)
            import                 :: dp
            character, intent(in)  :: norm
            integer  , intent(in)  :: m
            integer  , intent(in)  :: n
            real(dp) , intent(in)  :: a(lda, *)
            integer  , intent(in)  :: lda
            real(dp) , intent(out) :: work(*)
        endfunction dlange
    endinterface

    ! end BLAS-like ---------------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Auxiliary routines ----------------------------------------------------------------------------------

    ! end Auxiliary routines ------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------

endmodule lapack_interface