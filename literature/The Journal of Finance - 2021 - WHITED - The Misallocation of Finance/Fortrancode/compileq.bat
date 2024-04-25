erase *.exe
ifort /fast /O3 /c twoK_mod.f90
ifort /Qopenmp  compstat.f90 twoK_mod.f90 twoK_solve.f90 twoK_sim.f90 twoK_equil.f90 twoK_agg.f90 twoK_makeaggregates.f90 twoK_makemoments.f90 utilities.f90 /link /stack:8000000000
:: ifort /traceback /Qopenmp  compstat.f90 twoK_mod.f90 twoK_solve.f90 twoK_sim.f90 twoK_equil.f90 twoK_agg.f90 twoK_makeaggregates.f90 twoK_makemoments.f90 utilities.f90 /link /stack:8000000000
erase *.mod
erase *.obj
