module lapack_interface
    use iso_fortran_env, only: real32, real64
    implicit none

    integer, parameter, private :: sp = real32
    integer, parameter, private :: dp = real64

    interface
        subroutine dgeqp3(m  , n   , a    , lda , jpvt,     &
                          tau, work, lwork, info)
            import   :: dp
            integer  :: m
            integer  :: n
            real(dp) :: A(lda, *)
            integer  :: lda
            integer  :: jpvt(*)
            real(dp) :: tau(*)
            real(dp) :: work(*)
            integer  :: lwork
            integer  :: info
        endsubroutine dgeqp3
    endinterface

endmodule lapack_interface