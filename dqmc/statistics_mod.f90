module statistics_mod
    use numbertypes
    implicit none


    contains

        real(dp) function vector_avg(x, n)
            real(dp), intent(in) :: x(n)
            integer , intent(in) :: n

            vector_avg = sum(x) / n


        endfunction vector_avg


        subroutine jackknife(x, n, avg, err)
            real(dp), intent(in)  :: x(n)
            integer , intent(in)  :: n
            real(dp), intent(out) :: avg, err

            avg = vector_avg(x, n)
            err = sqrt(sum((x - avg) ** 2) / (n*(n-1)))

            
        endsubroutine jackknife


endmodule statistics_mod