# QM2Thermo
## Post-processing Thermoelectric Parameters from Quantum Mechanical calculations Output


## Components (Fortran 90 codes)
- generate_stencil.f90
- group_velocity.f90
- chemical_potential.f90
- Seebeck_analysis.f90


## Input files
- parameter.txt
- phononDOS.dat (option)
- lambda (option)
- a2F.dos* (option) from 
- Note: lambda and a2F.dos: QE format (Data from Abinit, etc. can be used if it matches the QE format.)


## WIEN2k output files (test: WIEN2k ver.12 and ver.16. LDA, PBE, WC, PBEsol or TB-mBJ)
- wien.dos1
- wien.energy
- wien.kgen
- wien.klist
- wien.struct
- Note: Code processing flow: WIEN2k (output files) -> generate_stencil.f90 (output: f*.dat) -> group_velocity.f90 (output: AKK.DATA) -> chemical_potential.f90 (output: apot.dat) -> Seebeck_analysis.f90 (output: Seebeck_analysis.dat)


## Preparation: Compilation
1. sudo apt update
2. sudo apt -y install gfortran build-essential make
3. make


## Usage
1. Calculated with symmetry 1 (P1 structure): After the calculation with SCF (either GGA-PBE or TB-mBJ), run the calculation with as many k-points as possible (e.g., Si 58 x 58 x 58 mesh) without shifting in DOS.
2. Extract the necessary files from WIEN2k and put them in the same file. (Optional: Calculate phonon DOS, lambda, a2F.dos* in QE and put them in the same file.)
3. Edit parameter.txt and enter the conditions you want to calculate.
4. The calculation is run with "bash ./run.sh".
5. (Optional: run plot_Seebeck.gpl to display the figure)

## Test
- Ubuntu 18.04 LTS or Later
- WSL, Windows 10 or Later
- Intel: Core-i3 or Later


## Reference
- [1] lambda and a2F.dos* data: https://github.com/nguyen-group/QE-SSP/tree/master/gr/alpha/reference)
- [2] quantum ESPRESSO tutorial: http://www.cmpt.phys.tohoku.ac.jp/~koretsune/SATL_qe_tutorial/elphon.html


## Citation
1. Journal version
```
@article{sato2024seebeck,
  author    = {Sato, Hiroshi and Miyazaki, Hiroshi and Nishino, Yoshikazu and others},
  title     = {Quantitative Evaluation of Seebeck Coefficient using Linearized Boltzmann Transport Equation for Fe₂VAl-Based Compounds},
  journal   = {Journal of Phase Equilibria and Diffusion},
  volume    = {45},
  pages     = {397--415},
  year      = {2024},
  doi       = {10.1007/s11669-024-01086-y},
  url       = {https://doi.org/10.1007/s11669-024-01086-y}
}
```
2. arXiv version
```

```
