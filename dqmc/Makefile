fc := gfortran
#fc_flags := -Wextra -Wall -pedantic -g -fbounds-check -pg
fc_flags := -Ofast
bmult := bmult_mod.o
#bmult := bmultexact_mod.o
mods := numbertypes.o simulate_mod.o equalgreens_mod.o customla_mod.o checkerboard_mod.o metropolisratios_mod.o $(bmult) simulationsetup_mod.o measurements_mod.o statistics_mod.o convenientla_mod.o

all: main.out

main.out: main.f90 numbertypes.o simulate_mod.o
	$(fc) $(fc_flags) main.f90 $(mods) -llapack -lblas -o main.out


convenientla_mod.o: convenientla_mod.f90 numbertypes.o customla_mod.o
	$(fc) $(fc_flags) -c convenientla_mod.f90

customla_mod.o: customla_mod.f90 numbertypes.o
	$(fc) $(fc_flags) -c customla_mod.f90

numbertypes.o: numbertypes.f90
	$(fc) $(fc_flags) -c numbertypes.f90

bmult_mod.o: bmult_mod.f90 numbertypes.o checkerboard_mod.o simulationsetup_mod.o customla_mod.o
	$(fc) $(fc_flags) -c bmult_mod.f90

checkerboard_mod.o: checkerboard_mod.f90 numbertypes.o
	$(fc) $(fc_flags) -c checkerboard_mod.f90

equalgreens_mod.o: equalgreens_mod.f90 numbertypes.o simulationsetup_mod.o $(bmult)
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

simulationsetup_mod.o: simulationsetup_mod.f90 numbertypes.o checkerboard_mod.o convenientla_mod.o
	$(fc) $(fc_flags) -c simulationsetup_mod.f90

statistics_mod.o: statistics_mod.f90 numbertypes.o
	$(fc) $(fc_flags) -c statistics_mod.f90

bmultexact_mod.o: bmultexact_mod.f90 numbertypes.o simulationsetup_mod.o customla_mod.o
	$(fc) $(fc_flags) -c bmultexact_mod.f90

clean:
	rm -rf *.o *.mod *.out
