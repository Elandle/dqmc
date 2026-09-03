module blas_lapack_interface
    use stduse
    implicit none

    ! -----------------------------------------------------------------------------------------------------
    ! Scalar operations -----------------------------------------------------------------------------------

    interface
        real(sp) function scabs1(z)
            import                  :: sp
            complex(sp), intent(in) :: z
        endfunction scabs1
    endinterface

    interface
        real(dp) function dcabs1(z)
            import                  :: dp
            complex(dp), intent(in) :: z
        endfunction dcabs1
    endinterface

    ! end Scalar operations -------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Level 1 BLAS: vector ops ----------------------------------------------------------------------------
    
    interface
        subroutine dcopy(n, dx, incx, dy, incy)
            import                 :: dp
            integer , intent(in)   :: n
            real(dp), intent(in)   :: dx(*)
            integer , intent(in)   :: incx
            real(dp), intent(out)  :: dy(*)
            integer , intent(in)   :: incy
        endsubroutine dcopy
    endinterface

    interface
        subroutine dscal(n, da, dx, incx)
            import                  :: dp
            integer , intent(in)    :: n
            real(dp), intent(in)    :: da
            real(dp), intent(inout) :: dx(*)
            integer , intent(in)    :: incx
        endsubroutine dscal
    endinterface

    interface
        subroutine dswap(n, dx, incx, dy, incy)
            import                  :: dp
            integer , intent(in)    :: n
            real(dp), intent(inout) :: dx(*)
            integer , intent(in)    :: incx
            real(dp), intent(inout) :: dy(*)
            integer , intent(in)    :: incy
        endsubroutine dswap
    endinterface

    interface
        subroutine daxpy(n   , da, dx, incx, dy,   &
                         incy)
            import                  :: dp
            integer , intent(in)    :: n
            real(dp), intent(in)    :: da
            real(dp), intent(in)    :: dx(*)
            integer , intent(in)    :: incx
            real(dp), intent(inout) :: dy(*)
            integer , intent(in)    :: incy
        endsubroutine daxpy
    endinterface

    ! end Level 1 BLAS: vector ops ------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Level 2 BLAS: matrix-vector ops ---------------------------------------------------------------------
    
    interface
        subroutine dger(m, n   , alpha, x  , incx,  &
                        y, incy, a    , lda)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            real(dp), intent(in)    :: alpha
            real(dp), intent(in)    :: x(*)
            integer , intent(in)    :: incx
            real(dp), intent(in)    :: y(*)
            integer , intent(in)    :: incy
            real(dp), intent(inout) :: a(lda, *)
            integer , intent(in)    :: lda
        endsubroutine dger
    endinterface

    ! end Level 2 BLAS: matrix-vector ops -----------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Level 3 BLAS: matrix-matrix ops ---------------------------------------------------------------------
    
    interface
        subroutine dgemm(transa, transb, m  , n, k  ,   &
                         alpha , a     , lda, b, ldb,   &
                         beta  , c     , ldc)
            import                   :: dp
            character, intent(in)    :: transa
            character, intent(in)    :: transb
            integer  , intent(in)    :: m
            integer  , intent(in)    :: n
            integer  , intent(in)    :: k
            real(dp) , intent(in)    :: alpha
            real(dp) , intent(in)    :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(in)    :: b(ldb, *)
            integer  , intent(in)    :: ldb
            real(dp) , intent(in)    :: beta
            real(dp) , intent(inout) :: c(ldc, *)
            integer  , intent(in)    :: ldc
        endsubroutine dgemm
    endinterface

    interface
        subroutine dtrmm(side, uplo , transa, diag, m,      &
                         n   , alpha, a     , lda , b,      &
                         ldb)
            import                   :: dp
            character, intent(in)    :: side
            character, intent(in)    :: uplo
            character, intent(in)    :: transa
            character, intent(in)    :: diag
            integer  , intent(in)    :: m
            integer  , intent(in)    :: n
            real(dp) , intent(in)    :: alpha
            real(dp) , intent(in)    :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(inout) :: b(ldb, *)
            integer  , intent(in)    :: ldb
        endsubroutine dtrmm
    endinterface

    ! end Level 3 BLAS: matrix-matrix ops -----------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


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

endmodule blas_lapack_interface