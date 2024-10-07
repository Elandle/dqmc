fc := gfortran
#fc_flags := -Wextra -Wall -pedantic -g -fbounds-check -pg
fc_flags := -Ofast

tests = customlatest.out checkerboardtest.out checkerboardtest.out

all: main.out $(tests)

main.out: main.f90 numbertypes.o simulate_mod.o
	$(fc) $(fc_flags) main.f90 numbertypes.o simulate_mod.o equalgreens_mod.o customla_mod.o asvqrd_mod.o checkerboard_mod.o metropolisratios_mod.o bmult_mod.o simulationsetup_mod.o measurements_mod.o statistics_mod.o -llapack -lblas -o main.out


customlatest.out: customlatest.f90 customlatests_mod.o customla_mod.o numbertypes.o
	$(fc) $(fc_flags) customlatest.f90 customlatests_mod.o customla_mod.o numbertypes.o convenientla_mod.o -llapack -lblas -o customlatest.out

checkerboardtest.out: checkerboardtest.f90 customla_mod.o checkerboard_mod.o checkerboardtests_mod.o
	$(fc) $(fc_flags) checkerboardtest.f90 numbertypes.o checkerboard_mod.o checkerboardtests_mod.o customla_mod.o convenientla_mod.o -llapack -lblas -o checkerboardtest.out

customlatests_mod.o: customlatests_mod.f90 customla_mod.o numbertypes.o convenientla_mod.o
	$(fc) $(fc_flags) -c customlatests_mod.f90

checkerboardtests_mod.o: checkerboardtests_mod.f90 numbertypes.o customla_mod.o checkerboard_mod.o convenientla_mod.o
	$(fc) $(fc_flags) -c checkerboardtests_mod.f90

convenientla_mod.o: convenientla_mod.f90 numbertypes.mod customla_mod.mod
	$(fc) $(fc_flags) -c convenientla_mod.f90

customla_mod.o: customla_mod.f90 numbertypes.o
	$(fc) $(fc_flags) -c customla_mod.f90

numbertypes.o: numbertypes.f90
	$(fc) $(fc_flags) -c numbertypes.f90

asvqrd_mod.o: asvqrd_mod.f90 customla_mod.o numbertypes.o
	$(fc) $(fc_flags) -c asvqrd_mod.f90

bmult_mod.o: bmult_mod.f90 numbertypes.o checkerboard_mod.o simulationsetup_mod.o customla_mod.o
	$(fc) $(fc_flags) -c bmult_mod.f90

checkerboard_mod.o: checkerboard_mod.f90 numbertypes.o
	$(fc) $(fc_flags) -c checkerboard_mod.f90

equalgreens_mod.o: equalgreens_mod.f90 numbertypes.o simulationsetup_mod.o asvqrd_mod.o bmult_mod.o
	$(fc) $(fc_flags) -c equalgreens_mod.f90

measurements_mod.o: measurements_mod.f90 numbertypes.o
	$(fc) $(fc_flags) -c measurements_mod.f90

metropolisratios_mod.o: metropolisratios_mod.f90 numbertypes.o simulationsetup_mod.o
	$(fc) $(fc_flags) -c metropolisratios_mod.f90

simulate_mod.o: simulate_mod.f90 numbertypes.o simulationsetup_mod.o equalgreens_mod.o metropolisratios_mod.o measurements_mod.o statistics_mod.o
	$(fc) $(fc_flags) -c simulate_mod.f90

measurements_mod.o: measurements_mod.f90 numbertypes.o simulationsetup_mod.o statistics_mod.o
	$(fc) $(fc_flags) -c measurements_mod.f90

statistics_mod.o: statistics_mod.f90 numbertypes.o
	$(fc) $(fc_flags) -c statistics_mod.f90

simulationsetup_mod.o: simulationsetup_mod.f90 numbertypes.o checkerboard_mod.o
	$(fc) $(fc_flags) -c simulationsetup_mod.f90

statistics_mod.o: statistics_mod.f90 numbertypes.o
	$(fc) $(fc_flags) -c statistics_mod.f90



clean:
	rm -rf *.o *.mod *.out
