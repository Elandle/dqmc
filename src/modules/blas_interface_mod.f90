module blas_interface
    use iso_fortran_env, only: real32, real64
    implicit none

    integer, parameter, private :: sp = real32
    integer, parameter, private :: dp = real64

    interface

        subroutine dgemm(transa, transb, m  , n, k  ,   &
                         alpha , a     , lda, b, ldb,   &
                         beta  , c     , ldc)
            import    :: dp
            character :: transa
            character :: transb
            integer   :: m
            integer   :: n
            integer   :: k
            real(dp)  :: alpha
            real(dp)  :: a(lda, *)
            integer   :: lda
            real(dp)  :: b(ldb, *)
            integer   :: ldb
            real(dp)  :: beta
            real(dp)  :: c(ldc, *)
            integer   :: ldc
        endsubroutine dgemm
    endinterface

endmodule blas_interface