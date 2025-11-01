module utilities
    use numbertypes
    implicit none

    contains

        integer function sgn(x)
            !
            ! Returns the sign of x as an integer.
            !
            !     1 if x >= 0
            !    -1 if x <  0
            !
            real(dp), intent(in) :: x

            if (x .ge. 0.0_dp) then
                sgn = 1
            else
                sgn = -1
            endif
        endfunction sgn

        integer function del(i, j)
            integer, intent(in) :: i
            integer, intent(in) :: j

            if (i .eq. j) then
                del = 1
            else
                del = 0
            endif
        endfunction del

endmodule utilities