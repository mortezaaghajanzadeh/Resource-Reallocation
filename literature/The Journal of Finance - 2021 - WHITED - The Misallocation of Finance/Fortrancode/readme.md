To generate the data for Table IA.VI in the internet appendix, you will need the following files:

Compilation file: 

Makefile

Program files:

datagenMPI.f90
twoK_mod.f90
twoK_agg.f90
twoK_equil.f90
twoK_makeaggregates.f90
twoK_makemoments.f90
twoK_momentgen.f90
twoK_sim.f90
twoK_solve.f90
utilities.f90

Instructions: The data generation program was run on the Comet supercomputer using MPI. The output then needs to be compiled using the Linux cat utility. The actual regressions are run using the following Matlab file: 

regress.m

To generate the figures in Figure IA.1, you will need the following files:

Compilation file:

compileq.bat

Program files:

compstat.f90 
twoK_mod.f90
twoK_agg.f90
twoK_equil.f90
twoK_makeaggregates.f90
twoK_makemoments.f90
twoK_momentgen.f90
twoK_sim.f90
twoK_solve.f90
utilities.f90

The baseline parameter values are in: 

estfil.txt
otherparam.txt




