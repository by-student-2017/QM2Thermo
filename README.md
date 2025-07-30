# LBT-TETRA (QM2Thermo)
## QM2Thermo: Post-processing Thermoelectric Parameters from Quantum Mechanical calculations Output
- Currently, only WIEN2k + QE(option data) is supported, but in the future we plan to make it compatible with various first-principles band calculations. We also plan to support GPUs with OpenACC in the future.
- Functions other than chemical potential and Seebeck coefficient, such as carrier concentration Nc, were created for future expansion and should be used only as a reference.
- Various function expansions will be fully implemented after new project applications are approved.
- H. Sato et al., J. Phase Equilib. Diffus. 45, 397–415 (2024).: [https://doi.org/10.1007/s11669-024-01086-y](https://doi.org/10.1007/s11669-024-01086-y)

## Preparation: Compilation
```
git clone https://github.com/by-student-2017/LBT-TETRA.git
cd LBT-TETRA
sudo apt update && sudo apt -y install gfortran build-essential make dos2unix
dos2unix *
make
```
Our group uses a mixture of operating systems, including Linux, Windows (WSL1 or WLS2), and Mac, so we use "dos2unix" to convert line breaks to Linux format.


## Recompilation
To recompile after changing the code, use:
```
make clean
make
```
Simply using "make" is fine, but if you are concerned, please use "make clean", which will also delete modules created with Fortran 90 code. (Note: Everything will be deleted, so do not use it in other directories.)


## Advanced Compilation
There are compile options in the Makefile. If you want to customize it, remove the # from "FLAGS" or "DEBUG" or add the necessary options. An example of the optimization option is:
```
FLAGS = -O2 -ftree-vectorize -ffp-contract=fast -fno-math-errno -march=native
```
The "-O3" option is not used here because it is an optimization that, in principle, has side effects. Therefore, "-O3" is not usually recommended, but it may be used in cases where the code is coded with side effects in mind, such as OpenMX (Professor Ozaki, the developer of OpenMX, is a genius, so he can do it. I, who am not very knowledgeable, could not do it). This code (LBT-TETRA) prioritizes readability as much as possible to make it easier for others to improve.  

Use DEBUG when you encounter problems improving your code. For example:
```
DEBUG = -g -Wall -fcheck=all -fcheck=bounds -fcheck=pointer -fbacktrace
```
As you can see from the previous examples, you only need to remove the "#". Of course, you can add or remove options.


## Usage
1. Calculated with symmetry 1 (P1 structure): After the calculation with SCF (either GGA-PBE or TB-mBJ), run the calculation with as many k-points as possible (e.g., Si Nk=34^3 (i.e., 58 x 58 x 58 k-points)) without shifting in DOS.
   - **Note**: When the lattice constants (a, b, c) are long or in the case of a supercell, the number of k-points can be reduced in proportion to "1/a x 1/b x 1/c". For example, in the FCC-type 2x2x2 supercell in the paper, Nk=34^3 (i.e., 34 x 34 x 34 k-points).)
2. Create a file named wien and put the WIEN2k output file into it. (Optional: In QE, for phonon DOS, lambda, and a2F.dos*, put it in the same file as the Fortran code *.f90.) You can use other file names by rewriting run.sh.
3. Edit parameter.txt and enter the conditions you want to calculate. (Note: ince this is Fortran, the read position is fixed. Please write the value within the provided range.)
4. The calculation is performed using the following command:
```
bash ./run.sh
```
5. To generate the plots using Gnuplot, run the following commands:
```
gnuplot plot_Seebeck.gpl
gnuplot plot_cp.gpl
gnuplot plot_dos.gpl
```
Graphs are also output in *.png format. If it doesn't work properly on Linux, you can install gnuplot on Windows and set up the environment, then just double-click the *.gpl file to make it work. (see *.png in "Ref_Si_TB-mBJ_dope")


### Visualization Scripts Overview
| Script Name         | Data File Used              | Visualization Target         | Description                                                                 |
|---------------------|-----------------------------|-------------------------------|-----------------------------------------------------------------------------|
| `plot_Seebeck.gpl`  | `Seebeck_analysis.dat`      | Seebeck Coefficient           | Visualizes the temperature dependence of the Seebeck coefficient related to thermoelectric effects. |
| `plot_cp.gpl`       | `apot.dat`                  | Chemical Potential            | Plots the variation of chemical potential in the material.                 |
| `plot_dos.gpl`      | `wien.dos1`                 | Density of States (DOS)       | Displays the electronic density of states for band structure analysis.     |

For DOS, you need to write the EF described in wien.dos1 into the script gnuplot plot_dos.gpl. The EF part is on line 22 of plot_dos.pgl shown below (written in [Ry] units). 
```
EF = 0.34363
```
As the second line of wien.dos1 states "EF = 0.34363", we can confirm that the same value is specified for EF on line 22 of plot_dos.pgl.   
- Ref_Si_TB-mBJ_dope shows an advanced use case where the Fermi level is shifted (using DEF in parameter.txt) to simultaneously plot two files: the conduction band (entering a plus value in DEF) and the valence band (entering a minus value in DEF). Note that "_plus" and "_minus" are manually added to the end of the output file (*.dat), and the gnuplot script (*_dope.gpl) also has more text. The VEC value is listed in Seebeck_analysis.dat (output file).

---

## Input files
first-principles codes, particularly for Seebeck coefficient and electron-phonon coupling analysis.
| File Name         | Purpose and Notes                                                                                      |
|------------------|--------------------------------------------------------------------------------------------------------|
| `parameter.txt`   | Used by `Seebeck_analysis.f90` to select calculation methods and input parameters. Also used in `chemical_potential.f90` for DEF-related settings. |
| `phononDOS.dat`   | *(Optional)* Contains phonon density of states (DOS) vs. energy (in eV). Can be sourced from any first-principles code. |
| `lambda`          | *(Optional)* Currently supported only in Quantum ESPRESSO (QE) format.                                |
| `a2F.dos*`        | *(Optional)* Currently supported only in Quantum ESPRESSO (QE) format. The first column (frequency/energy) and second column (a2Fdos_total) are essential; other columns are read but not used. The "*" in a2F.dos* corresponds to the number of data listed in the lambda file (an integer). |

lambda and a2F.dos: QE format (Data from Abinit, etc. can be used if it matches the QE format.)


## WIEN2k output files (test: WIEN2k ver.12 and ver.16. LDA, PBE, WC, PBEsol or TB-mBJ)
This document provides a summary of key files used in WIEN2k calculations and their respective purposes.
| File Name       | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `wien.dos1`     | Used to obtain data related to the Fermi level.                            |
| `wien.energy`   | Used to extract energy values at each k-point.                             |
| `wien.kgen`     | Used to retrieve the actual k-points utilized in the calculation.          |
| `wien.klist`    | Provides detailed information about the k-points (e.g., coordinates).      |
| `wien.struct`   | Contains lattice constants, angles, and symmetry data of the crystal.      |
In reality, fewer files would be sufficient, but this is done for convenience. If wien.klist is available, wien.kgen can actually be coded so that it is not necessary.


## Output files

| File Name       | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `AKK.DATA`      | Energy [eV], group velocity^2 [(m/s)^2], group velocity^2 * Total DOS [(m/s)^2 * (states/eV)], Total DOS [states/eV], Cumulative DOS [states/eV] are listed. |
| `apot.dat`      | The temperature [K] and chemical potential [eV] are listed. |
- AKK.DATA: AKK.DATA is output by group_velocity.f90, and is important data as it is used later by chemical_potential.f90 and Seebeck_analysis.f90. This information is essential when rewriting the code to calculate various physical quantities in Seebeck_analysis.f90.
- apot.dat: Since Seebeck_analysis.f90 uses this, the comments are placed at the end. If you want to change the list of temperatures to be calculated, you need to rewrite chemical_potential.f90 and Seebeck_analysis.f90 so that they correspond. This part is not dynamically allocated, so you need to match the number of data and change the array and information in the code.


## File Deletion Rules for Recalculation

If calculation result files already exist, the corresponding calculations will be skipped.  
To force recalculation, please delete the relevant files as described below.  
(You can find more details by looking at the run.sh Bash file.)

| Calculation Target     | Files to Delete         | Notes                                                                 |
|------------------------|-------------------------|-----------------------------------------------------------------------|
| Change in k-point mesh | `f*.dat`                | Example: `f001.dat`, `f002.dat`, etc. Delete files corresponding to the k-points. |
| Change in DEF          | `apot.dat`              | Delete this file if the DEF (structure details) has been modified.   |
| Other parameter changes| No deletion needed (`parameter.txt`) | Changes in `parameter.txt` alone will trigger recalculation automatically. |

> Note: Calculations will only run if the corresponding result files are not present.

---

## Test
- Ubuntu 18.04 LTS or Later
- WSL, Windows 10 or Later
- Intel: Core(TM)-i3-2100 or Later


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

## Components (Fortran 90 codes)
- generate_stencil.f90
- group_velocity.f90
- chemical_potential.f90
- Seebeck_analysis.f90


## Seebeck Coefficient Analysis Workflow

This document describes the workflow for computing the Seebeck coefficient using WIEN2k and a series of Fortran codes based on Allen's theory.

### Workflow Summary

| Code                   | Input File(s)                                                                               | Output File(s)                                            | Description                                           |
|------------------------|---------------------------------------------------------------------------------------------|-----------------------------------------------------------|-------------------------------------------------------|
| WIEN2k                 | -                                                                                           | wien.dos1, wien.energ, wien.kgen, wien.klist, wien.struct | Generate electronic structure data for the material.  |
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

## Troubleshooting
- The data reading position differs depending on the version of WIEN2k, so the Fortran 90 code needs to be rewritten.
- Calculations that have data files will be skipped, so if you want to calculate again, delete the file.

---

## How to Develop It Further

This section focuses on potential improvements to `Seebeck_analysis.f90`. Other parts of the codebase are already highly complete, especially for WIEN2k, and generally do not require rewriting.

### Key Areas for Improvement

- **Main Targets**:  
  The `seebeck_coefficient` and `get_tau` functions at the end of the main program are central.  
  Tracing `get_tau` reveals the part that calculates relaxation time using **Matthiessen's law**.  
  -> *Improving this calculation is highly recommended.*

- **Phonon Data Usage**:  
  Currently, only **Phonon DOS** is used.  
  -> *Further accuracy can be achieved by incorporating data from each q-point in the phonon spectrum.*

- **Alternative Method**:  
  As illustrated in the figure (not included here), another approach is to:
  - Use only the **chemical potential data** from `apot.dat`.
  - Adjust the Fermi level using **BoltzTraP** or similar tools.
  - Use only a subset of the data for analysis.

### Compatibility with Other First-Principles Codes

To make the code compatible with outputs from other DFT packages:

- **Preferred Format**:  
  Output data should follow the format of `wien.energy` and `wien.klist`.

- **Alternative Approach**:  
  Modify the data loading section in `group_velocity.f90`.  
  This module currently:
  - Retrieves energy values from `wien.energy`.
  - Loads k-point information from `wien.klist`, although this can also be extracted from `wien.energy`.

> Note: The current implementation uses `wien.energy` and `wien.klist` exactly as in the original publication, without modification.

### Developer Notes

- `generate_stencil.f90` (for HCP systems) has room for improvement.  
  However, due to the complexity and cost of development, modifications are **not recommended** unless you are an experienced developer. In practice, calculations involving element substitution, the introduction of various defects, and distortion result in the P1 symmetry, which is the highest, so there is little benefit to making it an HCP other than for research purposes.

---
