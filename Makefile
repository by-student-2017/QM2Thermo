# Fortran compiler settings
FC = gfortran
FLAGS = -O2 # -ftree-vectorize -ffp-contract=fast -fno-math-errno -march=native
DEBUG = 
#DEBUG = -g -Wall -fcheck=all -fcheck=bounds -fcheck=pointer -fbacktrace

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
