module printing_mod
    use stduse
    implicit none

    character(len=*), parameter :: dmatrixfmt = "(f17.8)"
    character(len=*), parameter :: dvectorfmt = "(f17.8)"
    character(len=*), parameter :: ivectorfmt = "(i6)"
    character(len=*), parameter :: imatrixfmt = "(i6)"

    interface print_matrix
        module procedure :: print_dmatrix
        module procedure :: print_imatrix
    endinterface

    interface print_vector
        module procedure :: print_dvector
        module procedure :: print_ivector
    endinterface

    
    contains

        subroutine print_dmatrix(A, ounit, message)
            real(dp)        , intent(in)           :: A(:, :)
            integer         , intent(in)           :: ounit
            character(len=*), intent(in), optional :: message
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            if (present(message)) then
                write (ounit, "(a)") message
            endif

            do i = 1, m
                do j = 1, n
                    write(unit=ounit, fmt=dmatrixfmt, advance="no") A(i, j)
                enddo
                write(unit=ounit, fmt="(a)") ""
            enddo
        endsubroutine print_dmatrix

        subroutine print_dvector(v, ounit, advance)
            real(dp), intent(in) :: v(:)
            integer , intent(in) :: ounit
            logical , optional   :: advance
            
            integer :: m, i
            logical :: adv

            if (present(advance)) then
                adv = advance
            else
                adv = .false.
            endif

            m = size(v)

            if (adv) then
                do i = 1, m
                    write(unit=ounit, fmt=dvectorfmt) v(i)
                enddo
            else
                do i = 1, m
                    write(unit=ounit, fmt=dvectorfmt, advance="no") v(i)
                enddo
                write(unit=ounit, fmt="(a)") ""
            endif
        endsubroutine print_dvector

        subroutine print_ivector(v, ounit, advance)
            integer, intent(in) :: v(:)
            integer, intent(in) :: ounit
            logical, optional   :: advance
            
            integer :: m, i
            logical :: adv

            if (present(advance)) then
                adv = advance
            else
                adv = .false.
            endif

            m = size(v)

            if (adv) then
                do i = 1, m
                    write(unit=ounit, fmt=ivectorfmt) v(i)
                enddo
            else
                do i = 1, m
                    write(unit=ounit, fmt=ivectorfmt, advance="no") v(i)
                enddo
                write(unit=ounit, fmt="(a)") ""
            endif
        endsubroutine print_ivector

        subroutine print_imatrix(A, ounit, message)
            integer         , intent(in)           :: A(:, :)
            integer         , intent(in)           :: ounit
            character(len=*), intent(in), optional :: message
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            if (present(message)) then
                write (ounit, "(a)") message
            endif

            do i = 1, m
                do j = 1, n
                    write(unit=ounit, fmt=imatrixfmt, advance="no") A(i, j)
                enddo
                write(ounit, "(a)") ""
            enddo
        endsubroutine print_imatrix
        
endmodule printing_mod