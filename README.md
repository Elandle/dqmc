# dqmc

This code (temporary name newQUEST) intends to be an updated version of QUEST (see https://github.com/Meromorphics/quest). While well performing, QUEST suffers from severe technical debt/feature creep/lack of documentation; leading it to be unmaintainable by new developers.
One of the primary goals of newQUEST is to be able to be maintained after years of nonmaintenance, even if that means sacrificing minor performance.

Currently newQUEST is in beta. This means it can actually run Monte Carlo (ie, properly sample Hubbard-Stratonavich fields, which means it has access to equal-time Green's functions and equal-time measurements can be easily taken).

newQUEST will be considered to be in version 1 once its single-thread (ie, no mpi or gpu) capabilities roughly match those of QUEST. Version 2 once mpi support has been integrated, and version 3 once gpu support has been integrated.
