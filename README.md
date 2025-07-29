## QM2Thermo: Post-processing Thermoelectric Parameters from Quantum Mechanical calculations Output


# Components (Fortran 90 codes)
- generate_stencil.f90
- group_velocity.f90
- chemical_potential.f90
- Seebeck_analysis.f90


# Input files
- parameter.txt
- phononDOS.dat (option)
- lambda (option) from QE format (Data from Abinit, etc. can be used if it matches the QE format.)
- a2F.dot* (option) from QE format (Data from Abinit, etc. can be used if it matches the QE format.)


# WIEN2k output files (test: WIEN2k ver.16. LDA, PBE, WC, PBEsol or TB-mBJ)
- wien.dos1
- wien.energy
- wien.kgen
- wien.klist
- wien.struct


# Preparation: Compilation
1. sudo apt update
2. sudo apt -y install gfortran build-essential make
3. make


# Usage
1. Calculated with symmetry 1 (P1 structure): After the calculation with SCF (either GGA-PBE or TB-mBJ), run the calculation with as many k-points as possible without shifting in DOS.
2. Extract the necessary files from WIEN2k and put them in the same file. (Optional: Calculate phonon DOS, lambda, a2F.dos* in QE and put them in the same file.)
3. Edit parameter.txt and enter the conditions you want to calculate.
4. The calculation is run with "bash ./run.sh".
5. (Optional: run plot_Seebeck.gpl to display the figure)

# Test
- Ubuntu 18.04 LTS or Later
- WSL, Windows 10 or Later
- Intel: Core-i3 or Later


# Reference
- [1] lambda and a2F.dos* data: https://github.com/nguyen-group/QE-SSP/tree/master/gr/alpha/reference)
- [2] quantum ESPRESSO tutorial: http://www.cmpt.phys.tohoku.ac.jp/~koretsune/SATL_qe_tutorial/elphon.html
