# LBT-TETRA (QM2Thermo)
## QM2Thermo: Post-processing Thermoelectric Parameters from Quantum Mechanical calculations Output
- Currently, only WIEN2k + "(option) EPW (from QE or Abinit, wannier90, etc) output" is supported, but in the future we plan to make it compatible with various first-principles band calculations. We also plan to support GPUs with OpenACC in the future.
- Functions other than chemical potential and Seebeck coefficient, such as carrier concentration Nc, were created for future expansion and should be used only as a reference.
- Various function expansions will be fully implemented after new project applications are approved.
- H. Sato et al., J. Phase Equilib. Diffus. 45, 397-415 (2024).: [https://doi.org/10.1007/s11669-024-01086-y](https://doi.org/10.1007/s11669-024-01086-y)

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
Simply using "make" is fine, but if you are concerned, please use "make clean", which will also delete modules created with Fortran 90 code. (Note: Do not use "make clean" in other directories as it will delete all *.exe and *.mod files in the directory.)


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
```
(On WIEN2k: (After SCF run) -> Tasks -> DOS -> x kgen)
Number of k-points: [0]
Shift k-mesh (if applicable) [No]
(For experts:...)
[58] [58] [58]
```
   - **Note**: When the lattice constants (a, b, c) are long or in the case of a supercell, the number of k-points can be reduced in proportion to "1/a x 1/b x 1/c". For example, in the FCC-type 2x2x2 supercell in the paper, Nk=34^3 (i.e., 34 x 34 x 34 k-points).)
2. Create a file named wien and put the WIEN2k output file into it. (Optional: In QE, for phonon DOS, lambda, and a2F.dos*, put it in the same file as the Fortran code *.f90.) You can use other file names by rewriting run.sh. If you use the "cd" and "ls" commands to look at the contents of the wien directory, you will see the following:
```
# inukai@PC-LAB3:/mnt/c/Users/inukai/sato_program/LBT-TETRA$ cd wien && ls && cd ../
cd wien && ls && cd ../
```
```
Si.dos1  Si.energy  Si.kgen  Si.klist  Si.struct
# Note: Although the names such as Si.dos1 have "Si" in them,
# they are renamed to wien by "run.sh" and placed in LBT-TETRA for use in the calculation.
```
   - **Option**: lambda, a2F.dos*, and phononDOS.dat should be placed in the LBT-TETRA directory. Only the necessary files are shown here. (Note: wien contains Si.dos1 Si.energy Si.kgen Si.klist Si.struct.)
```
# inukai@PC-LAB3:/mnt/c/Users/inukai/sato_program/LBT-TETRA$ ls
ls
```
```
a2F.dos1   a2F.dos2  a2F.dos4  a2F.dos6  a2F.dos8  lambda         wien
a2F.dos10  a2F.dos3  a2F.dos5  a2F.dos7  a2F.dos9  phononDOS.dat
```
3. Edit parameter.txt and enter the conditions you want to calculate. (Note: since this is Fortran, the read position is fixed. Please write the value within the provided range. Making sure to align the ":" and "!" will reduce errors.) For example:
```
!-----------------------:------------!--Memo---------------------------------------------------------------------------
DEF (Energy shift)  [eV]:  1.173000  ! This is good !
Base relaxation time [s]:  1.0       ! This is good !

!-----------------------:------------!--Memo---------------------------------------------------------------------------
DEF (Energy shift)  [eV]1.173000     ! This is bad ! (Note that it may be read as 0.173000.)
Base relaxation time [s]: 1.00e-14  !  This is bad ! ("!" is entered and an error occurs.)

!-----------------------:------------!--Memo---------------------------------------------------------------------------
DEF (Energy shift)  [eV]:        1.173000  ! This is bad ! (Note that it may be read as 1.17)
Base relaxation time [s] : 1.00e-14   ! This is bad ! (":" is entered and an error occurs.)

# As an aside, in the paper the basic relaxation time is calculated as 1.0.
# From this version onwards, the calculation of relaxation time is considered based on Alle's theory.
```
4. The calculation is performed using the following command:
```
bash ./run.sh
```
5. To generate the plots using Gnuplot, run the following commands:
```
gnuplot plot_Seebeck.gpl
gnuplot plot_cp.gpl
gnuplot plot_dos.gpl
gnuplot plot_ABGV2D.gpl
gnuplot plot_AB.gpl
```
Graphs are also output in *.png format. If it doesn't work properly on Linux, you can install gnuplot on Windows and set up the environment, then just double-click the *.gpl file to make it work. (see *.png in "Ref_Si_TB-mBJ_dope")


### Visualization Scripts Overview
| Script Name         | Data File Used              | Visualization Target         | Description                                                                 |
|---------------------|-----------------------------|-------------------------------|-----------------------------------------------------------------------------|
| `plot_Seebeck.gpl`  | `Seebeck_analysis.dat`      | Seebeck Coefficient           | Visualizes the temperature dependence of the Seebeck coefficient related to thermoelectric effects. |
| `plot_cp.gpl`       | `apot.dat`                  | Chemical Potential            | Plots the variation of chemical potential in the material.                 |
| `plot_dos.gpl`      | `wien.dos1`                 | Density of States (DOS)       | Displays the electronic density of states for band structure analysis.     |
| `plot_ABGV2D.gpl`   | `ABGV2D.dat`                | the spectrum A(E,T) and B(E,T), and <\|v\|^2 * DOS> | plot_ABGV2D.gpl visualizes the energy dependence of spectra A(E,T) and B(E,T) at the temperature specified by TEMP, and the electronic structure parameter <\|v\|^2 * DOS>. |
| `plot_AB.gpl`   | `Seebeck_analysis.dat`          | A(T) and B(T) | Plot the numerator A(T) and denominator B(T) as shown in Figure 12 in the paper. A(T) is the integral of A(E,T) with respect to the energy E, and B(T) is the integral of B(E,T) with respect to the energy E. |
- **plot_ABGV2D.gpl**: Please note that in this code, the relaxation time is explicitly specified and calculated, so the dimensions are multiplied by [s], unlike Figures 13-15 in the paper.
- **plot_AB.gpl**: Please note that in this code, the relaxation time is explicitly specified and calculated, so the dimensions are multiplied by [s], unlike Figures 12 in the paper.
- To specify the bands as shown in Figure 12-15, run the calculation in group_velocity.f90 using select_band_range.txt and output AKK.DATA containing only the specified bands.

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
| `lambda`          | *(Optional)* Currently supported only in EPW format.                                |
| `a2F.dos*`        | *(Optional)* Currently supported only in EPW format. The first column (frequency/energy) and second column (a2Fdos_total) are essential; other columns are read but not used. The "*" in a2F.dos* corresponds to the number of data listed in the lambda file (an integer). |
| `select_band_range.txt` | This is used to limit the range of the bands to be calculated, and is used in analyses such as Figures 13-15 in the paper. |


### WIEN2k output files (test: WIEN2k ver.12 and ver.16. LDA, PBE, WC, PBEsol or TB-mBJ)
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
If a calculation result file already exists, the corresponding calculation will be skipped. To force a recalculation, delete the corresponding file using the following procedure. (For details, see the run.sh Bash file.) You can run the calculation as usual with "bash ./run.sh".
| Calculation Target     | Files to Delete         | Notes                                                                 |
|------------------------|-------------------------|-----------------------------------------------------------------------|
| Change in k-point mesh | `cf*.dat`                | Delete only cfA1.dat file if change in k-point mesh or symmetry group. |
| Change in DEF          | `apot.dat`              | Delete apot.dat file if the DEF (structure details) has been modified.   |
| Other parameter changes| No delete files | Changing parameters other than DEF in `parameter.txt` does not require deleting the file. |
> Note: Calculations will only run if the corresponding result files are not present.

---

## Test
- Ubuntu 18.04 LTS or Later
- WSL, Windows 10 or Later
- Intel: Core(TM)-i3-2100 or Later


## References
- [1] lambda and a2F.dos* data: https://github.com/nguyen-group/QE-SSP/tree/master/gr/alpha/reference
- [2] quantum ESPRESSO tutorial: http://www.cmpt.phys.tohoku.ac.jp/~koretsune/SATL_qe_tutorial/elphon.html
- [3] T. Takeuchi, JTSJ (2011): https://doi.org/10.50972/thermoelectrics.8.1_17 (Open Access) (Japanese)
- [4] T. Takeuchi, JTSJ (2011): https://doi.org/10.50972/thermoelectrics.8.2_18 (Open Access) (Japanese)
- [5] T. Takeuchi, JTSJ (2012): https://doi.org/10.50972/thermoelectrics.8.3_27 (Open Access) (Japanese)
- [6] T. Takecuhi, JTSJ (2012): https://doi.org/10.50972/thermoelectrics.9.1_21 (Open Access) (Japanese)
- [7] T. Takeuchi, JTSJ (2012): https://doi.org/10.50972/thermoelectrics.9.2_21 (Open Access) (Japanese)
- [8] Y. Katsura, JTSJ (2014): https://doi.org/10.50972/thermoelectrics.10.3_20 (Open Access) (Japanese)
- [9] Y. Katsura, JTSJ (2014): https://doi.org/10.50972/thermoelectrics.11.1_18 (Open Access) (Japanese)
- [10] Y. Katsura, JTSJ (2014): https://doi.org/10.50972/thermoelectrics.11.2_19 (Open Access) (Japanese)
- [11] Empirical estimation of thermal conductivity: https://github.com/houzf/empirical_thermal_conductivity
- [12] Phonon density of states for Si: https://lampz.tugraz.at/~hadley/ss1/phonons/dos/si_phonon_dos.html
- [13] Thermal Conductivity of the Elements: https://srd.nist.gov/jpcrdreprint/1.3253100.pdf
- [14] Thermophysical Properties of Fluids Group: https://webbook.nist.gov/chemistry/fluid/
- [15] eXtremes of heat conduction: https://users.mrl.illinois.edu/cahill/mrs_sympx_f11.pdf
### Note
- If the bulk modulus and density in parameter.txt are > 0, the sound speed is calculated, and the thermal conductivity is calculated using the formula in "Empirical estimation of thermal conductivity".
- phononDOS.dat calculates the constant volume specific heat Cv_DOS and Debye temperature Theta_D. If the bulk modulus or density in parameter.txt are <= 0, it calculates the average sound speed from Cv_DOS and the thermal conductivity.
- The temperature dependence of the thermal conductivity is calculated using the Cezairliyan treatment of "Thermal Conductivity of Elements", where the code assumes that Tm = Debye temperature and km is the thermal conductivity at the Debye temperature.
- The constant A in the Slack model is derived from the formula commented in the python code in "Empirical estimation of thermal conductivity." This formula is closer to the thermal conductivity calculated from phonon DOS for Si than A = 3.1e-6, so we decided to use it.
- For the same bulk modulus, Poisson's ratio, and density, the Clarke model yields values approximately 20-40% lower than the Cahill model. This can be seen from the theoretical formula.

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
| generate_stencil.f90   | WIEN2k output files                                                                         | cf*.dat                                                    | Convert WIEN2k output into stencil data.              |
| group_velocity.f90     | cf*.dat, WIEN2k output files, select_band_range.txt                                         | AKK.DATA                                                  | Calculate group velocities from stencil data.         |
| chemical_potential.f90 | AKK.DATA, parameter.txt, WIEN2k output files                                                | apot.dat                                                  | Determine the chemical potential.                     |
| Seebeck_analysis.f90   | cf*.dat, apot.dat, AKK.DATA, parameter.txt, WIEN2k output files (optional: lambda, a2F.dos*) | Seebeck_analysis.dat, ABGV2D.dat                         | Compute the Seebeck coefficient using Allen's theory. A(E,T), B(E,T), <\|v\|^2 x DOS> are used to calculate the Seebeck coefficient using Allen's theory.|

### Code Descriptions

#### 1. WIEN2k
- **Purpose**: Generate electronic structure data for the material.
- **Output**: Band structure, density of states, and other relevant files.

#### 2. generate_stencil.f90
- **Purpose**: Convert WIEN2k output into stencil data.
- **Input**: WIEN2k output files.
- **Output**: `cf*.dat` files containing stencil information.

#### 3. group_velocity.f90
- **Purpose**: Calculate group velocities from stencil data.
- **Input**: `cf*.dat`, WIEN2k output files, and `select_band_range.txt`.
- **Output**: `AKK.DATA` containing velocity information.

#### 4. chemical_potential.f90
- **Purpose**: Determine the chemical potential.
- **Input**: `AKK.DATA`, `parameter.txt`, and WIEN2k output files.
- **Output**: `apot.dat` with chemical potential values.

#### 5. Seebeck_analysis.f90
- **Purpose**: Compute the Seebeck coefficient using Allen's theory.
- **Input**: `cf*.dat`, `apot.dat`, `AKK.DATA`, `parameter.txt`, and WIEN2k output files (optional: `lambda`, `a2F.dos*` in EPW format).
- **Output**: `Seebeck_analysis.dat` with final Seebeck coefficient results. `ABGV2D.dat` with the Spectrum A(E,T) and B(E,T), and Electronic structure parameter, <\|v\|^2 x DOS> results.

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

## Precision Considerations for `REAL(8)` (Double Precision)

### Overview

`REAL(8)` in Fortran corresponds to **IEEE 754 double precision**, which provides approximately **15–17 significant decimal digits**, typically around **16 digits**.

This level of precision is suitable for many scientific and engineering applications, but it has limitations that developers should be aware of when performing numerical computations.

---

### Machine Epsilon

- The **machine epsilon** for `REAL(8)` is approximately `2.220D-16`.
- This is the smallest value such that `1.0D0 + ε ≠ 1.0D0`.
- Differences smaller than this may be lost due to rounding errors.

---

### Practical Implications

When comparing floating-point numbers or checking for convergence, avoid direct equality checks. Instead, use a threshold to determine approximate equality:

```fortran
IF (ABS(x - y) < 1.0D-14) THEN
    ! x and y are considered approximately equal
END IF
```

---

# Allen-Type Seebeck Coefficient with Explicit Group Velocity

## Theoretical Background

The Seebeck coefficient $$\( S(T) \)$$ is calculated based on energy-dependent transport contributions, including group velocity, relaxation time, and electronic density of states. It is given by:

$$
S(T) = -\frac{1}{eT} \cdot 
\frac{
\int_{-\infty}^{\infty} (\varepsilon - \mu) \tau(\varepsilon) v^2(\varepsilon) g(\varepsilon) 
\left( -\frac{\partial f}{\partial \varepsilon} \right) d\varepsilon
}{
\int_{-\infty}^{\infty} \tau(\varepsilon) v^2(\varepsilon) g(\varepsilon) 
\left( -\frac{\partial f}{\partial \varepsilon} \right) d\varepsilon
}
$$

### Definitions

- $$\( S(T) \)$$: Seebeck coefficient (V/K)  
- $$\( T \)$$: Absolute temperature (K)  
- $$\( e \)$$: Elementary charge (C)  
- $$\( \varepsilon \)$$: Electron energy (eV)  
- $$\( \mu \)$$: Chemical potential (Fermi level)  
- $$\( f(\varepsilon) \)$$: Fermi-Dirac distribution  
- $$\( g(\varepsilon) \)$$: Electronic density of states  
- $$\( \tau(\varepsilon) \)$$: Energy-dependent relaxation time  
- $$\( v(\varepsilon) \)$$: Electron group velocity

## Relaxation Time from EPW

To evaluate $$\( \tau(\varepsilon) \)$$ from first-principles, the EPW (Electron-Phonon Wannier) code is used to calculate the Eliashberg spectral function $$\( \alpha^2F(\omega) \)$$. The relaxation time is computed using the following expression:

$$
\frac{1}{\tau(\varepsilon)} = 2\pi \int_0^\infty \alpha^2F(\omega) \left[ 1 + 2n(\omega) + f(\varepsilon - \omega) - f(\varepsilon + \omega) \right] d\omega
$$

Where:
- $$\( \alpha^2F(\omega) \)$$: Electron-phonon spectral function  
- $$\( n(\omega) \)$$: Bose–Einstein distribution  
- $$\( f(\varepsilon \pm \omega) \)$$: Fermi–Dirac distribution  

To activate this calculation within `epw.in`, include the directive:

```fortran
tau_vs_e = .true.
```

---

# Seebeck Coefficient Calculation with Kernel Filtering

## Overview

This document outlines the theoretical framework for computing the Seebeck coefficient $$\( S(T) \)$$, incorporating energy-dependent transport quantities and a generalized kernel (filter) function $$\( K(\varepsilon) \)$$ to selectively weight contributions across the energy spectrum.

---

## Theoretical Background

The Seebeck coefficient is defined by:

$$
S_K(T) = -\frac{1}{eT} \cdot 
\frac{
\int_{-\infty}^{\infty} (\varepsilon - \mu) \tau(\varepsilon) v^2(\varepsilon) g(\varepsilon) K(\varepsilon) \, d\varepsilon
}{
\int_{-\infty}^{\infty} \tau(\varepsilon) v^2(\varepsilon) g(\varepsilon) K(\varepsilon) \, d\varepsilon
}
$$

### Definitions

| Symbol                   | Description                             |
|--------------------------|-----------------------------------------|
| $$\( S(T) \)$$           | Seebeck coefficient (V/K)               |
| $$\( T \)$$              | Absolute temperature (K)                |
| $$\( e \)$$              | Elementary charge (C)                   |
| $$\( \varepsilon \)$$    | Electron energy (eV)                    |
| $$\( \mu \)$$            | Chemical potential (Fermi level)        |
| $$\( \tau(\varepsilon) \)$$ | Energy-dependent relaxation time     |
| $$\( v(\varepsilon) \)$$ | Electron group velocity                 |
| $$\( g(\varepsilon) \)$$ | Electronic density of states (DOS)      |
| $$\( K(\varepsilon) \)$$ | Kernel (filter) function                |
| $$\( f(\varepsilon) \)$$ | Fermi-Dirac distribution                |

---

## Kernel Function Design

The kernel function $$\( K(\varepsilon) \)$$ modulates the contribution of states to the transport integrals. Examples include:

| Type             | Formula                                              | Purpose                          |
|------------------|------------------------------------------------------|----------------------------------|
| Fermi window     | $$\( K(\varepsilon) = -\frac{\partial f}{\partial \varepsilon} \)$$ | Thermal broadening (standard)    |
| Gaussian         | $$\( K(\varepsilon) = \exp\left[ -\frac{(\varepsilon - \mu)^2}{2\sigma^2} \right] \)$$ | Emphasize local states           |
| Lorentzian       | $$\( K(\varepsilon) = \frac{1}{(\varepsilon - \mu)^2 + \gamma^2} \)$$ | Quasi-resonant states            |
| Moment-weighted  | $$\( K(\varepsilon) = (\varepsilon - \mu)^n \cdot \left( -\frac{\partial f}{\partial \varepsilon} \right) \)$$ | Higher-order transport moments   |

---

## Implementation Notes

- Numerical integration can be performed using quadrature schemes (e.g. Gauss-Legendre)
- Material-specific parameters for $$\( \tau \)$$, $$\( v \)$$, and $$\( g \)$$ may be obtained from DFT
- For custom kernels, normalization may be required for stability

---

## Use Cases

- Thermoelectric material screening  
- Quasi-local band contributions in complex oxides  
- Kernel optimization via machine learning  

---

## Group Velocity Definition

The group velocity $$\( v(\varepsilon) \)$$ is derived from the band structure and defined as:

$$
v(\varepsilon) = \frac{1}{\hbar} \cdot \frac{d\varepsilon}{dk}
$$

- $$\( \hbar \)$$: Reduced Planck constant  
- $$\( \frac{d\varepsilon}{dk} \)$$: Gradient of the energy band dispersion

---

## Energy-Dependent Relaxation Time via $$\( \alpha^2F(\omega) \)$$ (Implemented in Seebeck_analysis.f90)

The relaxation time $$\( \tau(\varepsilon) \)$$ can be calculated from the Eliashberg spectral function $$\( \alpha^2F(\omega) \)$$ using:

$$
\frac{1}{\tau(\varepsilon)} = 2\pi \int_0^{\infty} \alpha^2F(\omega) \left[
1 + 2n(\omega) + f(\varepsilon - \omega) - f(\varepsilon + \omega)
\right] d\omega
$$

Where:

- $$\( \alpha^2F(\omega) \)$$: Electron-phonon spectral function (from `a2F.dos`)  
- $$\( n(\omega) \)$$: Bose-Einstein distribution  
- $$\( f(\varepsilon \pm \omega) \)$$: Fermi-Dirac distribution  
- $$\( \omega \)$$: Phonon frequency

The file `a2F.dos` typically contains discretized $$\( \omega \) vs. \( \alpha^2F(\omega) \)$$ data, used in numerical evaluation of this integral.  
**Note:** This formulation is consistent with the implementation in [EPW](https://epw.gitlab.io/epw/), and the current code adopts the same theoretical approach to compute \( \tau(\varepsilon) \) based on electron–phonon interactions. The methodology allows direct compatibility with EPW output and facilitates integration with transport models such as the Allen-type Seebeck coefficient.

---

## Matthiessen’s Rule

If multiple scattering mechanisms contribute to relaxation, the total scattering rate is approximated via Matthiessen’s rule:

$$
\frac{1}{\tau_{\text{total}}(\varepsilon)} = 
\frac{1}{\tau_{\text{phonon}}(\varepsilon)} +
\frac{1}{\tau_{\text{impurity}}(\varepsilon)} +
\frac{1}{\tau_{\text{boundary}}(\varepsilon)} + \cdots
$$

This additive inverse relation allows individual mechanisms (phonon, impurity, grain boundaries, etc.) to be independently modeled and summed.

---

## Common Relaxation Time Approximations

| Scattering Type        | Approximation                          | Energy Dependence              |
|------------------------|----------------------------------------|--------------------------------|
| Acoustic phonon        | $$\( \tau(\varepsilon) \propto \varepsilon^{-1/2} \)$$ | Typical at high $$\( T \)$$        |
| Ionized impurity       | $$\( \tau(\varepsilon) \propto \varepsilon^{3/2} \)$$  | Relevant in doped semiconductors |
| Optical phonon         | Step-like or constant above threshold  | Depends on phonon energy       |
| Constant Approximation | $$\( \tau(\varepsilon) = \tau_0 \)$$       | Simplified CRTA (common fallback) |

## Relaxation Time Approximations for Transport Calculations

| Scattering Type             | Approximate Expression                                                     | Energy Dependence           | Temperature Dependence     | Applicability             |
|----------------------------|-----------------------------------------------------------------------------|-----------------------------|-----------------------------|---------------------------|
| **Acoustic phonon**        | $$\( \tau(\varepsilon) \propto \varepsilon^{-1/2} \)$$                         | Weak                        | Mild                        | Metals, high-T            |
| **Ionized impurity**       | $$\( \tau(\varepsilon) \propto \|\varepsilon - \mu\|^n \)$$                      | Strong (user-defined $$\( n \)$$) | Weak or none               | Doped semiconductors      |
| **Constant (CRTA)**        | $$\( \tau(\varepsilon) = \tau_0 \)$$                                           | None                        | None                        | Simple models (use 1.0D-14 [s])  |
| **Eliashberg-based**       | See: $$\( \alpha^2F(\omega) \)$$ integration                                   | Fully resolved              | Fully resolved              | First-principles accurate |
| **$$\( \lambda \)$$-based Allen approx.** | $$\( \tau(T) \approx \frac{1}{\pi \lambda k_B T} \)$$                     | None (Fermi-level only)     | Linear in $$\( T \)$$          | Metals                    |
| **DOS-based (electronic)** | $$\( \tau(\varepsilon) \propto \frac{1}{g(\varepsilon) T} \)$$                | Moderate                    | Linear in $$\( T \)$$          | Generic materials         |
| **Phonon scattering**      | $$\( \tau(\varepsilon) \propto \frac{1}{\|\varepsilon - \mu\| \cdot T} \)$$      | Linear near band edges      | Linear                      | Semiconductors, general   |
| **Phonon DOS-based**       | $$\( \tau^{-1} \propto \int D_{\text{ph}}(\omega) \cdot \frac{1 + n(\omega,T)}{\omega} \, d\omega \)$$ | Indirect via phonons        | Strong via $$\( n(\omega) \)$$ | Full phonon spectrum      |

> Notes:
> - $$\( g(\varepsilon) \)$$: electronic density of states  
> - $$\( D_{\text{ph}}(\omega) \)$$: phonon density of states  
> - $$\( n(\omega, T) \)$$: Bose-Einstein distribution

---

### Implementation Tips

- Use analytical forms or interpolated data to model DOS or phonon spectra.
- Normalize relaxation time expressions for dimensionless comparison if needed.
- Numerical integration required for phonon DOS-based model.

> Note: Accurate modeling requires evaluating $$\( \tau(\varepsilon) \)$$ either from first-principles $$\( \alpha^2F(\omega) \)$$, or using fitted functional forms consistent with experimental data.

---

## Relaxation Time Estimation Using Coupling Constant $$\( \lambda \)$$

When the full Eliashberg spectral function $$\( \alpha^2F(\omega) \)$$ is unavailable, the energy-averaged **electron-phonon coupling constant** $$\( \lambda \)$$ can be used to estimate scattering rates via simplified models.

The coupling constant is defined as:

$$
\lambda = 2 \int_0^{\infty} \frac{\alpha^2F(\omega)}{\omega} d\omega
$$

Once $$\( \lambda \)$$ is known (often extracted from experiments or DFT calculations), the **electron-phonon scattering rate** $$\( \tau^{-1} \)$$ can be approximated using:

### Simplified Allen Formula:

$$
\frac{1}{\tau(\varepsilon)} \approx \pi \lambda k_B T
$$

This model assumes:

- Quasiparticles near the Fermi level $$\( \varepsilon \approx \mu \)$$  
- Weak energy dependence in the scattering  
- $$\( T \)$$-linear behavior typical in metals

It provides a temperature-dependent relaxation time:

$$
\tau(T) \approx \frac{1}{\pi \lambda k_B T}
$$

### Notes:

- This approximation is most valid at **high temperatures** (above Debye temperature) and for **isotropic metals**.
- It cannot resolve energy-resolved transport and should be avoided if detailed band-dependent scattering rates are needed.
- Useful for quick estimates or initial screening of thermoelectric candidates.

---

## Coupling Constant vs. Full $$\( \alpha^2F(\omega) \)$$

| Method                        | Input Requirement           | Resolution       | Energy Dependence | Accuracy |
|------------------------------|-----------------------------|------------------|-------------------|----------|
| Full Eliashberg Function     | `a2F.dos` from DFPT         | High (fine grid) | Yes               | ⭐⭐⭐⭐     |
| $$\( \lambda \)$$-based Estimate | Single scalar $$\( \lambda \)$$ | Low      | No                | ⭐⭐       |

---

# Mode-Resolved Electron-Phonon Scattering with Kernel Filtering

## Overview

This module implements a refined phonon-mode filtering and scattering kernel designed to improve the accuracy of relaxation time calculations in electron-phonon coupling analysis. It includes:

- **Mode-resolved phonon weighting filter** $$\( w_{\text{mode}}(\omega) \)$$
- **Kernel-based energy-frequency matching logic** $$\( K(\varepsilon, \omega) \)$$
- Support for Eliashberg function construction and selection rules

---

## Components

### 1. Phonon Mode Filter: `w_mode`

Each phonon mode is assigned a weight $$\( w_{\text{mode}}(\omega_i) \)$$ based on its frequency category:

| Frequency Region (ω in eV)     | Mode Type      | Weight $$\( w_{\text{mode}} \)$$ |
|-------------------------------------|----------------|-------------------------------|
| $$\( \omega < 0.136058 \)$$             | Acoustic       | 1.0                           |
| $$\( 0.01 \leq \omega < 0.408174 \)$$   | Intermediate   | 0.5                           |
| $$\( \omega \geq 0.408174 \)$$          | Optical        | 0.1                           |

These weights are initialized in the `read_phonon_dos` routine and refined in `generate_w_mode()`.

> **Note**: Additional logic for more granular filtering is available and can be toggled for precise control over modal contributions.

---

### 2. Kernel Function: `kernel_selection(E, ω)`

This logical kernel enforces an energy–frequency matching rule in scattering calculations:

$$
K(\varepsilon, \omega) = 
\begin{cases}
1.0 & \text{if } \omega < |\varepsilon| \\
0.0 & \text{otherwise}
\end{cases}
$$

**Purpose**: Ensures that only phonon modes with energy below the electronic excitation $$\( \varepsilon \)$$ participate in scattering.

---

## Relaxation Time Evaluation

The relaxation time $$\( \tau(\varepsilon) \)$$ is computed as:

$$
\frac{1}{\tau(\varepsilon)} = \sum_{i = 1}^{N_{\text{ph}}} g_{\text{eff}}(\omega_i) \cdot w_{\text{mode}}(\omega_i) \cdot K(\varepsilon, \omega_i) \cdot F(\omega_i)
$$

Where:

- $$\( g_{\text{eff}}(\omega_i) \)$$: Mode-resolved coupling strength (default: 1.0)
- $$\( F(\omega_i) \)$$: Spectral term from phonon DOS or Eliashberg function
- $$\( K(\varepsilon, \omega_i) \)$$: Kernel selection result (1.0 or 0.0)

---

## Reference
- B. Xu, M. Di Gennaro, M. J. Verstraete, *Phys. Rev. B* **102**, 155128 (2020).: [https://doi.org/10.1103/PhysRevB.102.155128](https://doi.org/10.1103/PhysRevB.102.155128)

---

