# dqmc

This is my determinant quantum Monte Carlo (DQMC) code.
It was built with the intention of handling asymmetric hoppings between sites.
I am updating this code a bit, but I intend on making this code "legacy" when I finish writing a new (and object-oriented Fortran!) version.

# Compiling

Compilation is done using `cmake`.
I have got a `fpm.toml` from an older iteration if you prefer that, but I have no guarantee that it still works (I know it does not have any optimization flags).

## Basic Installation, Compilation, and Running
First obtain the code and `cd` into it:
```
    git clone https://github.com/Elandle/dqmc.git
    cd dqmc
```
For a first compilation, do a basic "native" build and compilation:
```
    cmake -B build
    cmake --build build
```
Now the main executable `./build/src/executables/main` is compiled.
You can move it to wherever you want to run simulations at.
Here I will make a folder called `simulation` and copy it over:
```
    mkdir simulation
    cp ./build/src/executables/main ./simulation
```
Now to get a simulation rolling.
Three files for input are required (an example is in the directory `example`): a checkerboard file, a bipartite file, and an input file.
(The bipartite file is a required "unrequired" file. Right now I am doing bipartite runs so I have not bothered yet to unrequire it. Just create a txt with `i 1` all the way down for `i = 1, 2, ..., N` to fake putting all sites in the same sublattice and ignore measurements requiring this).
An input file has the following structure:
```
    N           = <number of sites in lattice>
    L           = <number of imaginary time slices>
    nstab       = <maximum number of Green's function wraps to tune up to before recalculation>
    north       = <number of B matrices multiplied together before a QRP factorization in a new Green's function calculation>
    nbin        = <number of bins to put measurements into. Must divide nmeassweep (otherwise you will get wrong results which I have not implemented a check for yet)>
    nmeassweep  = <number of measurements to take>
    nskip       = <number of sweeps to skip between measurements>
    nequil      = <number of warmup sweeps>
    dtau        = <imaginary time step>
    U           = <electron interaction strength>
    mu          = <chemical potential>
    ckbfilename = <name of file to read the checkerboard from>
    bipfilename = <name of file to read bipartite from>
    outfilename = <name of file to output results into>
    debfilename = <name of file for debugging output>
    summary     = <line to print at the top of results file>
```
You can use as many spaces on a line as you want and comment the rest of a line using `#`.
Checkerboard and bipartite files are generated using Python code in the directory `ckb` (I have not yet implemented a graph datastructure and colouring algorithm in Fortran).
There are some examples in that folder, but the basic thing you do is create a matrix `K` of hoppings (`K[i, j] = tij = hopping site j to site i`), then call:
```
    check = ckb.ckb(K)
    check.saveckb(ckbname)
    check.savebipartite(bipname)
```
to save the checkerboard input file in the string held by `ckbname` and similarly for the bipartite input file for `bipname`.


Once you get the checkerboard and bipartite input files setup, put them in `simulation`, type up an input file, then run like:
```
    ./main <inputfile>
```
The code will announce a few things about the simulation (eg, how many imaginary time slices) then start up warmup and Green's function `nstab` tuning, and measuring.

## Optimization Options
I have a few `cmake` optimization options that affect how the code is compiled. They are:
```
    cmake -B build -Doptimization=nativehigh
    cmake -B build -Doptimization=puma
    cmake -B build -Doptimization=portable
    cmake -B build -Doptimization=debug
```
`nativehigh` is the default optimization option (what happens if you do not specify any). It applies many aggresive optimizations made for your cpu (eg, `-march=native`). Use this if you plan on running the code on your current computer. `debug` is slow and for debugging.
`puma` and `portable` are two options I wrote to run the code on computers I have access to. Some are part of a cluster called `puma` consisting of somewhat older intel cpus. Others are quite old (intel i3 or i5 3rd or 5th generation I believe) so I use `portable` for these. This is so I can compile on my main computer statically (eg, I do not need BLAS/LAPACK on the other computers) and send out the code to other computers to run.

## Note
It appears some LAPACK distributions don't have the routines `dlarscl2` and `dlascl2` (eg, my Mac has it but my Linux pc does not). If you get compilation errors due to this, go to `customla_mod.f90` and the `left_diaginvmult` and `left_diagmult` subroutines. Comment out those problematic calls and uncomment out the alternatives (just below).
