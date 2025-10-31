# dqmc

This code (temporary name newQUEST) intends to be an updated version of QUEST (see https://github.com/Meromorphics/quest). While well performing, QUEST suffers from severe technical debt/feature creep/lack of documentation; leading it to be unmaintainable by new developers.
One of the primary goals of newQUEST is to be able to be maintained after years of nonmaintenance, even if that means sacrificing minor performance.

Currently newQUEST is in beta. This means it can actually run Monte Carlo (ie, properly sample Hubbard-Stratonavich fields, which means it has access to equal-time Green's functions and equal-time measurements can be easily taken).

Here is a list of planned features for implementation:
* Openmp parallelization (including parallel checkerboard multiplication, rolling feeder Green's function calculation, unequal time Green's function calculation)
* MPI parallelization (multiple simulations at once)
* GPU implementation
* Enhanced geometry specification (through Python, with visualization of lattice)
* Unequal time Green's function calculations
* More measurements (catch up to around QUEST level)
* Delay update or submatrix update (?, maybe unnecessary complication)

# Compiling

newQUEST has moved to Fortran Package Manager (fpm https://fpm.fortran-lang.org/) for building from Make. This is still new to me, so I'll add instructions once I figure it out more. For now, if you have fpm installed, just type fpm build and the code will automagically compile. fpm makes a bunch of random files and folders, in which you can find a main file. This is the main thing to run. The plan is to move to CMake.

Note: it appears some LAPACK distributions don't have the routines dlarscl2 and dlascl2 (eg, my Mac has it but my Linux pc doesn't). If you get compilation errors due to this, go to customla_mod.f90 and the left_diaginvmult and left_diagmult subroutines. Comment out those problematic calls and uncomment out the alternatives (just below).
