program main
    use stduse
    use simulate_mod
    use iso_fortran_env, only: output_unit
    implicit none
    !
    ! Executable file that runs a simulation.
    ! At the moment this is being used as a simulation tester,
    ! so it is heavily under construction. The eventual plan is to
    ! move over to a bunch of similar example files that run
    ! simulations like this one.
    !
    ! Run like: ./main input.txt
    !
    type(Simulation)              :: S
    integer                       :: iunit, ilen
    character(len=:), allocatable :: ifile

    call get_command_argument(1, length=ilen)
    allocate(character(len=ilen) :: ifile)
    call get_command_argument(1, ifile)
    call setup_simulation_input(S, trim(ifile), iunit, output_unit)

    ! call setup_simulation_input(S, "input.txt", iunit, output_unit)
    S%h = 1
    ! Open debugging output file
    open(newunit=S%dunit, file=S%debfilename, action="write", status="replace")
    call simulate(S)

endprogram main