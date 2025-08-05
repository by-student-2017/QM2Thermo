# Fortran compiler settings
FC = gfortran
NCPUs := $(shell expr $(shell nproc) / 2)
FLAGS = -O2 -march=native -fopenmp -ftree-parallelize-loops=$(NCPUs) # -ftree-vectorize -ffp-contract=fast -fno-math-errno -floop-parallelize-all
DEBUG = 
#DEBUG = -g -Wall -fcheck=all -fcheck=bounds -fcheck=pointer -fbacktrace
##GPU: FC = nvfortran, FLAGS = -acc -ta=multicore -Minfo=acc

# Executables and their corresponding source files
all: group_velocity.exe chemical_potential.exe Seebeck_analysis.exe generate_stencil.exe

group_velocity.exe: group_velocity.f90
	$(FC) $(FLAGS) $(DEBUG) $< -o $@

chemical_potential.exe: chemical_potential.f90
	$(FC) $(FLAGS) $(DEBUG) $< -o $@

Seebeck_analysis.exe: Seebeck_analysis.f90
	$(FC) $(FLAGS) $(DEBUG) $< -o $@

generate_stencil.exe: generate_stencil.f90
	$(FC) $(FLAGS) $(DEBUG) $< -o $@

# Clean-up command
clean:
	rm -f *.exe *.mod
