module readinputfile_mod
    use numbertypes
    use iso_fortran_env, only: iostat_end
    implicit none
    !
    ! Basic usage:
    !
    ! First declare:
    ! type(parameter_values) :: param_values
    ! type(parameter_set)    :: param_set
    !
    ! To read in:
    ! readinputfile(fname, funit, ounit, param_values, param_set)
    !
    ! Where:
    !       fname is a string containing the input file name
    !       funit is an integer that will hold the file unit number of fname (fname will be opened then closed during the call)
    !       ounit is the file unit (an integer) where to print warnings (eg, ounit = output_unit from iso_fortran_env for the terminal)
    !       
    ! To print the contents:
    ! call printparams(param_values, ounit)
    !
    type parameter_values
        integer                       :: N
        integer                       :: L
        integer                       :: nstab
        integer                       :: north
        integer                       :: nbin
        integer                       :: nmeassweep
        integer                       :: nskip
        integer                       :: nequil
        real(dp)                      :: dtau
        real(dp)                      :: U
        real(dp)                      :: mu
        character(len=:), allocatable :: ckbfilename
    endtype parameter_values


    type parameter_set
        logical :: N           = .false.    
        logical :: L           = .false.
        logical :: nstab       = .false.
        logical :: north       = .false.
        logical :: nbin        = .false.
        logical :: nmeassweep  = .false.
        logical :: nskip       = .false.
        logical :: nequil      = .false.
        logical :: dtau        = .false.
        logical :: U           = .false.
        logical :: mu          = .false.
        logical :: ckbfilename = .false.
    endtype parameter_set


    contains


        subroutine readln(funit, iostat, line, maxlen)
            integer         ,              intent(in)              :: funit
            integer                      , intent(out)             :: iostat
            character(len=:), allocatable, intent(inout)           :: line
            integer         ,              intent(in)   , optional :: maxlen

            ! Cannot edit maxlen since it is optional, so using another integer to hold the value used
            integer :: actual_maxlen

            if (present(maxlen)) then
                actual_maxlen = maxlen
            else
                ! default maximum length of 100
                actual_maxlen = 100
            endif

            ! allocate if not already allocated
            if (.not. allocated(line)) then
                allocate(character(actual_maxlen) :: line)
            endif

            ! force line to be maxlen long
            if (len(line) .ne. actual_maxlen) then
                deallocate(line)
                allocate(character(actual_maxlen) :: line)
            endif

            read(unit=funit, fmt="(a)", iostat=iostat) line
            line = ltrim(line)
            

        endsubroutine readln


        subroutine rmcmt(line, cmt)
            character(len=:), allocatable, intent(inout)           :: line
            character(len=*),              intent(in)   , optional :: cmt

            character(len=:), allocatable :: actual_cmt
            integer :: loc

            if (present(cmt)) then
                actual_cmt = cmt
            else
                ! default comment of #
                actual_cmt = "#"
            endif

            loc = index(line, actual_cmt, back=.false.)

            if (loc .gt. 0) then
                line = ltrim(line(1 : loc-1))
            endif


            deallocate(actual_cmt)

        endsubroutine rmcmt


        logical function isasmt(line, asmtop)
            character(len=*), intent(in)           :: line
            character(len=*), intent(in), optional :: asmtop

            character(len=:), allocatable :: actual_asmtop

            if (present(asmtop)) then
                actual_asmtop = asmtop
            else
                ! default assignment operator of =
                actual_asmtop = "="
            endif

            if (index(line, actual_asmtop) .gt. 0) then
                isasmt = .true.
            else
                isasmt = .false.
            endif


        endfunction isasmt


        subroutine getasmnts(line, left, right, asmtop)
            character(len=*), intent(in)               :: line
            character(len=:), allocatable, intent(out) :: left
            character(len=:), allocatable, intent(out) :: right
            character(len=*), intent(in) , optional    :: asmtop

            character(len=:), allocatable :: actual_asmtop
            integer                       :: loc

            if (present(asmtop)) then
                actual_asmtop = asmtop
            else
                ! default assignment operator of =
                actual_asmtop = "="
            endif

            loc = index(line, actual_asmtop, back=.false.)

            if (loc .gt. 0) then
                left  = ltrim(line(1 : loc-1))
                right = ltrim(line(loc + len(actual_asmtop) :))
            else
                left  = ""
                right = ""
            endif


        endsubroutine getasmnts


        function ltrim(line)
            character(len=*), intent(in) :: line
            character(len=len(trim(adjustl(line)))) :: ltrim

            ltrim = trim(adjustl(line))

        endfunction ltrim


        subroutine asnparam(left, right, param_values, param_set)
            character(len=*)             , intent(in)    :: left
            character(len=*)             , intent(in)    :: right
            type(parameter_values)      , intent(inout) :: param_values
            type(parameter_set)         , intent(inout) :: param_set

            select case(left)
                case("N")
                    read(unit=right, fmt="(i20)") param_values%N
                    param_set%N = .true.
                case("L")
                    read(unit=right, fmt="(i20)") param_values%L
                    param_set%L = .true.
                case("nstab")
                    read(unit=right, fmt="(i20)") param_values%nstab
                    param_set%nstab = .true.
                case("north")
                    read(unit=right, fmt="(i20)") param_values%north
                    param_set%north = .true.
                case("nbin")
                    read(unit=right, fmt="(i20)") param_values%nbin
                    param_set%nbin = .true.
                case("nmeassweep")
                    read(unit=right, fmt="(i20)") param_values%nmeassweep
                    param_set%nmeassweep = .true.
                case("nskip")
                    read(unit=right, fmt="(i20)") param_values%nskip
                    param_set%nskip = .true.
                case("nequil")
                    read(unit=right, fmt="(i20)") param_values%nequil
                    param_set%nequil = .true.
                case("dtau")
                    read(unit=right, fmt="(f20.10)") param_values%dtau
                    param_set%dtau = .true.
                case("U")
                    read(unit=right, fmt="(f20.10)") param_values%U
                    param_set%U = .true.
                case("mu")
                    read(unit=right, fmt="(f20.10)") param_values%mu
                    param_set%mu = .true.
                case("ckbfilename")
                    param_values%ckbfilename = right
                    param_set%ckbfilename = .true.
            endselect


        endsubroutine asnparam


        subroutine printparams(param_values, ounit)
            type(parameter_values), intent(in) :: param_values
            integer               , intent(in) :: ounit

            write(unit=ounit, fmt="(a20, i20)")    "N = "          , param_values%N
            write(unit=ounit, fmt="(a20, i20)")    "L = "          , param_values%L
            write(unit=ounit, fmt="(a20, i20)")    "nstab = "      , param_values%nstab
            write(unit=ounit, fmt="(a20, i20)")    "north = "      , param_values%north
            write(unit=ounit, fmt="(a20, i20)")    "nbin = "       , param_values%nbin
            write(unit=ounit, fmt="(a20, i20)")    "nmeassweep = " , param_values%nmeassweep
            write(unit=ounit, fmt="(a20, i20)")    "nskip = "      , param_values%nskip
            write(unit=ounit, fmt="(a20, i20)")    "nequil = "     , param_values%nequil
            write(unit=ounit, fmt="(a20, f20.10)") "dtau = "       , param_values%dtau
            write(unit=ounit, fmt="(a20, f20.10)") "U = "          , param_values%U
            write(unit=ounit, fmt="(a20, f20.10)") "mu = "         , param_values%mu
            write(unit=ounit, fmt="(a20, a20)")    "ckbfilename = ", param_values%ckbfilename

        endsubroutine printparams


        subroutine asndefaults(param_values, param_set, ounit)
            type(parameter_values), intent(inout) :: param_values
            type(parameter_set)   , intent(inout) :: param_set
            integer                , intent(in)    :: ounit
            
            integer  :: L_default
            integer  :: nstab_default
            integer  :: north_default
            integer  :: nbin_default
            integer  :: nmeassweep_default
            integer  :: nskip_default
            integer  :: nequil_default
            real(dp) :: dtau_default
            real(dp) :: U_default
            real(dp) :: mu_default

            L_default          = 60
            nstab_default      = 5
            north_default      = 5
            nbin_default       = 32
            nmeassweep_default = 1000 * nbin_default
            nskip_default      = 5
            nequil_default     = 1000
            dtau_default       = 0.125
            U_default          = 0.0_dp
            mu_default         = 0.0_dp

            if (.not. param_set%N) then
                write(unit=ounit, fmt="(a)") "Warning: N not set. Simulation cannot run unless assigned."
            endif
            if (.not. param_set%L) then
                write(unit=ounit, fmt="(a)") "Warning: L not set. Setting to default value."
                param_values%L = L_default
                param_set%L    = .true.
            endif
            if (.not. param_set%nstab) then
                write(unit=ounit, fmt="(a)") "Warning: nstab not set. Setting to default value."
                param_values%nstab = nstab_default
                param_set%nstab    = .true.
            endif
            if (.not. param_set%north) then
                write(unit=ounit, fmt="(a)") "Warning: north not set. Setting to default value."
                param_values%north = north_default
                param_set%north    = .true.
            endif
            if (.not. param_set%nbin) then
                write(unit=ounit, fmt="(a)") "Warning: nbin not set. Setting to default value."
                param_values%nbin = nbin_default
                param_set%nbin    = .true.
            endif
            if (.not. param_set%nmeassweep) then
                write(unit=ounit, fmt="(a)") "Warning: nmeassweep not set. Setting to default value."
                param_values%nmeassweep = nmeassweep_default
                param_set%nmeassweep    = .true.
            endif
            if (.not. param_set%nskip) then
                write(unit=ounit, fmt="(a)") "Warning: nskip not set. Setting to default value."
                param_values%nskip = nskip_default
                param_set%nskip    = .true.
            endif
            if (.not. param_set%nequil) then
                write(unit=ounit, fmt="(a)") "Warning: nequil not set. Setting to default value."
                param_values%nequil = nequil_default
                param_set%nequil    = .true.
            endif
            if (.not. param_set%dtau) then
                write(unit=ounit, fmt="(a)") "Warning: dtau not set. Setting to default value."
                param_values%dtau = dtau_default
                param_set%dtau    = .true.
            endif
            if (.not. param_set%U) then
                write(unit=ounit, fmt="(a)") "Warning: U not set. Setting to default value."
                param_values%U = U_default
                param_set%U    = .true.
            endif
            if (.not. param_set%mu) then
                write(unit=ounit, fmt="(a)") "Warning: mu not set. Setting to default value."
                param_values%mu = mu_default
                param_set%mu    = .true.
            endif
            if (.not. param_set%mu) then
                write(unit=ounit, fmt="(a)") "Warning: ckbfilename not set. Simulation cannot run unless assigned."
            endif


        endsubroutine asndefaults


        subroutine readinputfile(fname, funit, ounit, param_values, param_set)
            character(len=*)      , intent(in)    :: fname
            integer               , intent(inout) :: funit
            integer               , intent(in)    :: ounit
            type(parameter_values), intent(inout) :: param_values
            type(parameter_set)   , intent(inout) :: param_set

            character(len=:), allocatable :: line
            character(len=:), allocatable :: left, right
            integer                       :: iostat

            open(newunit=funit, file=fname)

            ! read line-by-line until the end
            do while (iostat .ne. iostat_end)
                call readln(funit, iostat, line) ! read in line
                call rmcmt(line)                 ! get rid of comment (if present)
                ! not at the end of the file and the line isn't blank
                if (iostat .ne. iostat_end .and. len(line) > 0) then
                    ! what to do if this line is an assignment line
                    if (isasmt(line)) then
                        call getasmnts(line, left, right)                    ! read the left and right hand sides of assignment
                        call asnparam(left, right, param_values, param_set)  ! assign the left to be the right
                    endif
                endif
            enddo

            close(unit=funit)
    
            ! assign default values (if needed)
            call asndefaults(param_values, param_set, ounit)

        endsubroutine readinputfile



endmodule readinputfile_mod