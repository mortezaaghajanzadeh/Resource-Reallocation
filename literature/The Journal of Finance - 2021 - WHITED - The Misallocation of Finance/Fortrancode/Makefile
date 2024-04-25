FC = mpif90
FCFLAGS   = -g -O3 -xHOST -qopenmp -mkl -mcmodel medium -shared-intel
twoK: datagenMPI.f90 twoK_equil.o twoK_solve.o twoK_sim.o twoK_makemoments.o twoK_makeaggregates.o twoK_agg.o utilities.o  twoK_mod.o 
	$(FC) $(FCFLAGS) -o gen datagenMPI.f90 twoK_equil.o twoK_solve.o twoK_sim.o twoK_makemoments.o twoK_makeaggregates.o twoK_agg.o utilities.o  twoK_mod.o  
twoK_equil.o: twoK_equil.f90 twoK_solve.o twoK_sim.o twoK_makemoments.o twoK_makeaggregates.o  twoK_agg.o utilities.o  twoK_mod.o
	$(FC) $(FCFLAGS) -c  twoK_equil.f90 twoK_solve.o twoK_sim.o twoK_makemoments.f90 twoK_makeaggregates.o twoK_agg.o utilities.o twoK_mod.o
twoK_solve.o: twoK_solve.f90 twoK_sim.o twoK_makemoments.o utilities.o  twoK_mod.o
	$(FC) $(FCFLAGS) -c  twoK_solve.f90 twoK_sim.o twoK_makemoments.f90 utilities.o twoK_mod.o
twoK_sim.o: twoK_sim.f90 utilities.o  twoK_mod.o
	$(FC) $(FCFLAGS) -c  twoK_sim.f90 utilities.o twoK_mod.o
twoK_agg.o: twoK_agg.f90 utilities.o  twoK_mod.o
	$(FC) $(FCFLAGS) -c  twoK_agg.f90 utilities.o twoK_mod.o
twoK_makeaggregates.o: twoK_makeaggregates.f90 utilities.o  twoK_mod.o
	$(FC) $(FCFLAGS) -c  twoK_makeaggregates.f90 utilities.o twoK_mod.o
twoK_makemoments.o: twoK_makemoments.f90 utilities.o  twoK_mod.o 
	$(FC) $(FCFLAGS) -c  twoK_makemoments.f90 utilities.o twoK_mod.o 
utilities.o: utilities.f90 twoK_mod.o
	$(FC) $(FCFLAGS) -c  utilities.f90 twoK_mod.o
twoK_mod.o: twoK_mod.f90
	$(FC) $(FCFLAGS) -c twoK_mod.f90
