# dqmc

This code (temporary name newQUEST) intends to be an updated version of QUEST (see https://github.com/Meromorphics/quest). While well performing, QUEST suffers from severe technical debt/feature creep/lack of documentation; leading it to be unmaintainable by new developers.
One of the primary goals of newQUEST is to be able to be maintained after years of nonmaintenance, even if that means sacrificing minor performance.

Currently newQUEST is in beta. This means it can actually run Monte Carlo (ie, properly sample Hubbard-Stratonavich fields, which means it has access to equal-time Green's functions and equal-time measurements can be easily taken).

newQUEST will be considered to be in version 1 once its single-thread (ie, no mpi or gpu) capabilities roughly match those of QUEST. Version 2 once mpi support has been integrated, and version 3 once gpu support has been integrated.

# Compiling

newQUEST has moved to Fortran Package Manager (fpm https://fpm.fortran-lang.org/) for building from Make. This is still new to me, so I'll add instructions once I figure it out more. For now, if you have fpm installed, just type fpm build and the code will automagically compile. fpm makes a bunch of random files and folders, in which you can find a main file. This is the main thing to run.

Note: it appears some LAPACK distributions don't have the routines dlarscl2 and dlascl2 (eg, my Mac has it but my Linux pc doesn't). If you get compilation errors due to this, go to customla_mod.f90 and the left_diaginvmult and left_diagmult subroutines. Comment out those problematic calls and uncomment out the alternatives (just below).
