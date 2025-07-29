# LBT-TETRA (QM2Thermo)
## QM2Thermo: Post-processing Thermoelectric Parameters from Quantum Mechanical calculations Output
- Currently, only WIEN2k+QE is supported, but in the future we plan to make it compatible with various first-principles band calculations.
- H. Sato et al., J. Phase Equilib. Diffus. 45, 397–415 (2024). [https://doi.org/10.1007/s11669-024-01086-y](https://doi.org/10.1007/s11669-024-01086-y)


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


## References
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
---

## Seebeck Coefficient Analysis Workflow

This document describes the workflow for computing the Seebeck coefficient using WIEN2k and a series of Fortran codes based on Allen's theory.

### Workflow Summary

| Code                   | Input File(s)                                                                               | Output File(s)                                            | Description                                           |
|------------------------|---------------------------------------------------------------------------------------------|-----------------------------------------------------------|-------------------------------------------------------|
| WIEN2k                 | -                                                                                           | Band structure, DOS, and other electronic structure files | Generate electronic structure data for the material.  |
| generate_stencil.f90   | WIEN2k output files                                                                         | f*.dat                                                    | Convert WIEN2k output into stencil data.              |
| group_velocity.f90     | f*.dat, WIEN2k output files                                                                 | AKK.DATA                                                  | Calculate group velocities from stencil data.         |
| chemical_potential.f90 | AKK.DATA, parameter.txt, WIEN2k output files                                                | apot.dat                                                  | Determine the chemical potential.                     |
| Seebeck_analysis.f90   | f*.dat, apot.dat, AKK.DATA, parameter.txt, WIEN2k output files (optional: lambda, a2F.dos*) | Seebeck_analysis.dat                                      | Compute the Seebeck coefficient using Allen's theory. |

### Code Descriptions

#### 1. WIEN2k
- **Purpose**: Generate electronic structure data for the material.
- **Output**: Band structure, density of states, and other relevant files.

#### 2. generate_stencil.f90
- **Purpose**: Convert WIEN2k output into stencil data.
- **Input**: WIEN2k output files.
- **Output**: `f*.dat` files containing stencil information.

#### 3. group_velocity.f90
- **Purpose**: Calculate group velocities from stencil data.
- **Input**: `f*.dat` and WIEN2k output files.
- **Output**: `AKK.DATA` containing velocity information.

#### 4. chemical_potential.f90
- **Purpose**: Determine the chemical potential.
- **Input**: `AKK.DATA`, `parameter.txt`, and WIEN2k output files.
- **Output**: `apot.dat` with chemical potential values.

#### 5. Seebeck_analysis.f90
- **Purpose**: Compute the Seebeck coefficient using Allen's theory.
- **Input**: `f*.dat`, `apot.dat`, `AKK.DATA`, `parameter.txt`, and WIEN2k output files (optional: `lambda`, `a2F.dos*` in QE format).
- **Output**: `Seebeck_analysis.dat` with final Seebeck coefficient results.

---
