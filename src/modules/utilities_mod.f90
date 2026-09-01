module utilities
    use stduse
    implicit none

    abstract interface
        integer function seed_function(seed, i)
            integer, intent(in) :: seed
            integer, intent(in) :: i
        endfunction seed_function
    endinterface

    contains

        pure elemental integer function sgn(x)
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

        pure elemental integer function del(i, j)
            integer, intent(in) :: i
            integer, intent(in) :: j

            if (i .eq. j) then
                del = 1
            else
                del = 0
            endif
        endfunction del


        !
        ! Sets seed to a "random" (at least according to Fortran's intrinsic
        ! random_init(repeatable=.false.)) single integer seed.
        !
        subroutine rand_seed(seed, sfunction)
            integer                 , intent(out) :: seed
            procedure(seed_function), optional    :: sfunction

            real(dp) :: r

            call random_init(repeatable=.false., image_distinct=.true.)
            call random_number(r)
            seed = int((2.0_dp*r - 1.0_dp) * huge(seed))
            call set_seed(seed, sfunction)
        end subroutine rand_seed

        !
        ! Initializes the random number generator using seed.
        ! There is an optional sfunction argument with signature:
        !
        !       sfunction(seed, i)
        !
        ! where seed and i are intent(in) integer and returns an integer.
        ! This function bridges the gap between single integer seeds
        ! and Fortran's seed array (it populates the seed array in a loop
        ! over this, ie Fortran's seed(i) = sfunction(seed, i))).
        ! A "random" one made by a bunch of bit operations is provided by default
        ! and does not overflow.
        !
        subroutine set_seed(seed, sfunction)
            integer                 , intent(in) :: seed
            procedure(seed_function), optional   :: sfunction

            integer              :: n
            integer              :: i
            integer, allocatable :: put(:)

            call random_seed(size=n)
            allocate(put(n))

            if (present(sfunction)) then
                do i = 1, n
                    put(i) = sfunction(seed, i)
                end do
            else
                do i = 1, n
                    put(i) = seed_function_dflt(seed, i)
                end do
            end if
            call random_seed(put=put)
            deallocate(put)
        endsubroutine set_seed

        function seed_function_dflt(seed, i) result(j)
            integer, intent(in) :: seed
            integer, intent(in) :: i

            integer :: j

            ! Bunch of random bit stuff (avoids overflow)
            j = ieor(seed, ishftc(i, 7))
            j = ieor(j, shiftl(j, 5))
            j = ieor(j, shiftr(j, 17))
            j = ishftc(j, 9)
            j = ieor(j, not(shiftr(j, 11)))
            j = ieor(j, shiftl(j, 3))
            j = ishftc(j, -7)
            j = ieor(j, ishftc(i, 13))
            j = ieor(j, shiftr(j, 19))
            j = ieor(j, shiftl(j, 6))
            j = ishftc(j, 11)
            j = ieor(j, shiftr(j, 8))
            j = ieor(j, shiftl(j, 15))
        endfunction seed_function_dflt

endmodule utilities