program main
    use numbertypes
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
    type(Simulation) :: S
    integer          :: iunit
    call setup_simulation_input(S, "input.txt", iunit, output_unit)
    S%h = 1
    ! Open debugging output file
    open(newunit=S%dunit, file=S%debfilename, action="write", status="replace")
    call simulate(S)

endprogram main