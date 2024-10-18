module statistics_mod
    use numbertypes
    implicit none
    !
    ! Contains procedures for doing statistics on input data
    !
    contains

    
        real(dp) function vector_avg(x, n)
            !
            ! Returns the average of the entries of the length n vector x
            !
            real(dp), intent(in) :: x(n)
            integer , intent(in) :: n

            vector_avg = sum(x) / n


        endfunction vector_avg


        subroutine jackknife(x, n, avg, err)
            !
            ! Performs the jackknife method on the entries of the length n vector x.
            !
            ! The average obtained by the jackknife method will be stored in avg,
            ! and the error in err.
            !
            real(dp), intent(in)  :: x(n)
            integer , intent(in)  :: n
            real(dp), intent(out) :: avg, err

            avg = vector_avg(x, n)
            err = sqrt(sum((x - avg) ** 2) / (n*(n-1)))

            
        endsubroutine jackknife


endmodule statistics_mod