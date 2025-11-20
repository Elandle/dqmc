!> \brief Contains procedures and datatypes of reading in input files.
!!
!! While hardcoded for DQMC, the code can easily be adapted for other uses.
!!
!! Basic usage.
!!
!! First declare: <br>
!! `type(parameter_values) :: param_values` <br>
!! `type(parameter_set)    :: param_set`    <br>
!!
!! To read in values: <br>
!! `readinputfile(fname, funit, ounit, param_values, param_set)` <br>
!!
!! To print the contents: <br>
!! `call printparams(param_values, ounit)`
module readinputfile_mod
    use stduse
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

    !> \brief Holds read in parameter values for DQMC.
    type parameter_values
        integer                       :: N           !< Number of sites \f$N\f$ in lattice.
        integer                       :: L           !< Number of imaginary time slices \f$L\f$.
        integer                       :: nstab       !< Every `nstab`th time a time slice is advanced, the Green's functions \f$G_\sigma\f$ are computed from new as opposed to wrapping.
        integer                       :: north       !< Number of \f$B_\sigma\f$ matrices to multiply together before doing a \f$QRP\f$ factorisation in constructing a Green's function from scratch.
        integer                       :: nbin        !< How many bins to separate measurements into.
        integer                       :: nmeassweep  !< Total number of sweeps to perform a measurement on.
        integer                       :: nskip       !< How many sweeps to skip between measurements.
        integer                       :: nequil      !< How many warmup sweeps to perform.
        real(dp)                      :: dtau        !< Trotterization parameter \f$\Delta\tau\f$.
        real(dp)                      :: U           !< Interaction strength \f$U\f$.
        real(dp)                      :: mu          !< Chemical potential \f$\mu\f$.
        character(len=:), allocatable :: ckbfilename !< File name to read checkerboard from.
        character(len=:), allocatable :: outfilename !< File name to print results to.
        character(len=:), allocatable :: debfilename !< File name to print debugging information (if in debug mode) into.
    endtype parameter_values

    !> \brief Tracks whether or not a value has set for each parameter.
    !!
    !! Used when assigning default parameter values or stopping code execution if
    !! a vital parameter is not set. When a parameter is set (either by being read in
    !! or set to a default value) the corresponding `logical` variable in this type is
    !! set `.true.`.
    type parameter_set
        logical :: N           = .false.   !< `.true.` if `N` is set. `.false.` if not.
        logical :: L           = .false.   !< `.true.` if `L` is set. `.false.` if not.
        logical :: nstab       = .false.   !< `.true.` if `nstab` is set. `.false.` if not.
        logical :: north       = .false.   !< `.true.` if `north` is set. `.false.` if not.
        logical :: nbin        = .false.   !< `.true.` if `nbin` is set. `.false.` if not.
        logical :: nmeassweep  = .false.   !< `.true.` if `nmeassweep` is set. `.false.` if not.
        logical :: nskip       = .false.   !< `.true.` if `nskip` is set. `.false.` if not.
        logical :: nequil      = .false.   !< `.true.` if `nequil` is set. `.false.` if not.
        logical :: dtau        = .false.   !< `.true.` if `dtau` is set. `.false.` if not.
        logical :: U           = .false.   !< `.true.` if `U` is set. `.false.` if not.
        logical :: mu          = .false.   !< `.true.` if `mu` is set. `.false.` if not.
        logical :: ckbfilename = .false.   !< `.true.` if `ckbfilename` is set. `.false.` if not.
        logical :: outfilename = .false.   !< `.true.` if `outfilename` is set. `.false.` if not.
        logical :: debfilename = .false.   !< `.true.` if `debfilename` is set. `.false.` if not.
    endtype parameter_set

    contains

        !> \brief Reads the next line from a file.
        !!
        !! Reminder: in Fortran, files are streams.
        !! So if a line is read from a file (eg, by using the `read` intrinsic)
        !! that line is consumed in the stream: the next time a line is read,
        !! it is actually the next line in the file (to go back, an intrinsic
        !! such as `backspace` or `rewind` may be used).
        !!
        !! Reads the next line from file with unit `funit`, storing it in `line`
        !! with length `maxlen`. Leading spaces are removed, and `line` is padded
        !! with trailing spaces to fill up to `maxlen` length if necessary.
        !!
        !! Example: if the line being read is: </p>
        !! `nstab = 4` </p>
        !! Then if `readln` is called with `maxlen = 12`, the read in line is `"nstab = 4    "`. 
        !!
        !! \param[in]    funit  (`integer`)                      Unit of file to read from.
        !! \param[out]   iostat (`integer`)                      Status of the `read` statement on the file read from.
        !! \param[inout] line   (`character(len=:), allocatable`) String to allocate (`line` must be allocatable) read line into.
        !! \param[in]    maxlen (`integer, optional`)            Maximum length of a line to read and the exit length of `line`. Lines with data to be read in longer than `maxlen` being read into may result in errors. Default value: `100`.
        subroutine readln(funit, iostat, line, maxlen)
            integer                      , intent(in)              :: funit
            integer                      , intent(out)             :: iostat
            character(len=:), allocatable, intent(inout)           :: line
            integer                      , intent(in)   , optional :: maxlen

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

        !> \brief Removes a comment from a line (if a comment is present).
        !!
        !! Removes from `line` (must be `allocatable`) all text after and including `cmt` present in `line` and all leading and trailing spaces.
        !!
        !! For example, if `line = "  N = 12 # 12 site lattice"` and `cmt = "#"` (the default value),
        !! then after calling `rmcmt`,  `line = "N = 12"`.
        !!
        !! \param[inout] line (`character(len=:), allocatable`) String to remove comment and leading and trailing spaces from.
        !! \param[in]    cmt  (`character(len=:), optional`)    String denoting the start of a comment in `line`. Default value: `"#"`.
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
            else
                line = ltrim(line)
            endif

            deallocate(actual_cmt)
        endsubroutine rmcmt

        !> \brief Determines whether or not a line corresponds to an assignment.
        !!
        !! Determines whether or not `line` corresponds to an assignment by
        !! determining whether or not it has `asmtop` in it.
        !!
        !! For example, if `line = "N = 12"` and `asmtop = "="` (default),
        !! then `isasmt` would return `.true.`.
        !!
        !! \param[in] line    (`character(len=*)`)           String to determine whether or not an assignment is present in.
        !! \param[in] asmtop  (`character(len=*), optional`) Assignment operator to search for in string to determine whether it is an assignment or not. Default value: `"="`.
        !! \result    isasmt  (`logical`)                    `.true.` if `line` corresponds to an assignment (contains the string `asmtop`). `.false.` if not.
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

        !> \brief Removes leading and trailing blank spaces from an input string.
        !!
        !! Example: if `line = "  abc de   "`, then `ltrim(line)` returns
        !! `"abc de"`.
        !!
        !! \param[in] line   (`character(len=*)`)       String to remove leading and trailing blank spaced from.
        !! \result    ltrim  (`character(len=trimlen)`) `line` with leading and trailing blank spaces removed, with a length `trimlen` of `line` minus the amount of leading and trailing blank spaces.
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
                    param_set   %ckbfilename = .true.
                case("outfilename")
                    param_values%outfilename = right
                    param_set   %outfilename = .true.
                case("debfilename")
                    param_values%debfilename = right
                    param_set   %debfilename = .true.
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
            write(unit=ounit, fmt="(a20, a20)")    "outfilename = ", param_values%outfilename
            write(unit=ounit, fmt="(a20, a20)")    "debfilename = ", param_values%debfilename
        endsubroutine printparams

        subroutine asndefaults(param_values, param_set, ounit)
            type(parameter_values), intent(inout) :: param_values
            type(parameter_set)   , intent(inout) :: param_set
            integer               , intent(in)    :: ounit
            
            integer           :: L_default
            integer           :: nstab_default
            integer           :: north_default
            integer           :: nbin_default
            integer           :: nmeassweep_default
            integer           :: nskip_default
            integer           :: nequil_default
            real(dp)          :: dtau_default
            real(dp)          :: U_default
            real(dp)          :: mu_default
            character(len=10) :: outfilename_default

            L_default           = 60
            nstab_default       = 5
            north_default       = 5
            nbin_default        = 32
            nmeassweep_default  = 1000 * nbin_default
            nskip_default       = 5
            nequil_default      = 1000
            dtau_default        = 0.125
            U_default           = 0.0_dp
            mu_default          = 0.0_dp
            outfilename_default = "output.txt"

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
            if (.not. param_set%ckbfilename) then
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