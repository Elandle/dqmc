module blas_interface
    use stduse
    implicit none

    ! -----------------------------------------------------------------------------------------------------
    ! Scalar operations ----------------------------------------------------------------------------

    interface
        real(sp) function scabs1(z)
            import :: sp
            complex(sp), intent(in) :: z
        endfunction scabs1
    endinterface

    interface
        real(dp) function dcabs1(z)
            import                  :: dp
            complex(dp), intent(in) :: z
        endfunction dcabs1
    endinterface

    ! end Scalar operations ------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Level 1 BLAS: vector ops ----------------------------------------------------------------------------
    
    interface
        subroutine dcopy(n, dx, incx, dy, incy)
            import :: dp
            integer , intent(in)  :: n
            real(dp), intent(in)  :: dx(*)
            integer , intent(in)  :: incx
            real(dp), intent(out) :: dy(*)
            integer , intent(in)  :: incy
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
            real(dp), intent(inout) :: dy
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
            import   :: dp
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
            import    :: dp
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

endmodule blas_interface