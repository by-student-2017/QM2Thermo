# QE v6.8 + Elastic2020

## QE v6.8 Installation
1. sudo apt update
2. sudo apt -y install gfortran g++ build-essential make libopenblas-dev libopenmpi-dev libfftw3-dev
3. sudo apt -y install gnuplot ghostscript
3. wget https://github.com/QEF/q-e/archive/refs/tags/qe-6.8.tar.gz
4. tar xvf qe-6.8.tar.gz
5. cd q-e-qe-6.8
6. ./configure
7. make pwall
8. sudo make install

## Elastic2020 Installation
1. sudo apt update
2. sudo apt -y install python3 python3-numpy python3-matplotlib
3. cd $HOME
4. git clone https://github.com/by-student-2017/Elastic2020.git
5. cd Elastic2020
6. tar xvf adon_v1_0.tar.gz
7. cd SpaceGroups
8. make
9. cd ..
10. cp SpaceGroups/sgroup ./sgroup
11. echo 'export PATH="$HOME/Elastic2020:$PATH"' >> ~/.bashrc

## Calculation method (energy base and 2nd order)
1. run_elastic_energy_2nd.sh

## Calculation method (stress base and 2nd order)
1. run_elastic_stress_2nd.sh

## Calculation method (stress base and 3rd order) (Some confirmation required: need to check RII to N structure)
1. run_elastic_stress_3rd.sh

## Calculation method (energy base and 3rd order) (All need to be checked)
1. run_elastic_energy_3rd.sh

## Whether or not to display the diagram
- The figure is stored in a directory named Plot. If you want to see it immediately (*.sh), do the following:
```
# Choose whether to display the diagram: yes or no
show_plot_flag="yes"
```

## Delete old files
- You can delete some of the files using the command below. If you don't need any directories (folders), delete them too.
```
bash clean.sh
```

## input file (case.scf.in)
- Please use the format of the input file (case.scf.in) that is the base for QE. This is mainly because CELL_PARAMETERS will be changed and the location where it is written will also change.

## Test
- OS: Ubuntu 22.04 LTS (WLS2)
- Python: 3.10.12
- Numpy: 1.21.5
- Matplotlib: 3.5.1
- Scipy: 3.5.1

## Results: Energy and 2nd order: Diamond Si (Time: about 2 [min])
```
ElaStic_2nd.out found in Energy-vs-Strain. Displaying contents:
    The output of ElaStic code
    Today is Sat Aug 16 00:49:48 2025

    Symmetry of the second-order elastic constant matrix in Voigt notation.
    for, space group-number between 207 and 230, Cubic I structure.

               C11     C12     C12      0       0       0
               C12     C11     C12      0       0       0
               C12     C12     C11      0       0       0
                0       0       0      C44      0       0
                0       0       0       0      C44      0
                0       0       0       0       0      C44

    Elastic constant (stiffness) matrix in GPa:


      180.9        43.5        43.5         0.0         0.0         0.0
       43.5       180.9        43.5         0.0         0.0         0.0
       43.5        43.5       180.9         0.0         0.0         0.0
        0.0         0.0         0.0        68.5         0.0         0.0
        0.0         0.0         0.0         0.0        68.5         0.0
        0.0         0.0         0.0         0.0         0.0        68.5

    Elastic compliance matrix in 1/GPa:


    0.00610    -0.00118    -0.00118     0.00000     0.00000     0.00000
   -0.00118     0.00610    -0.00118     0.00000     0.00000     0.00000
   -0.00118    -0.00118     0.00610     0.00000     0.00000     0.00000
    0.00000     0.00000     0.00000     0.01459     0.00000     0.00000
    0.00000     0.00000     0.00000     0.00000     0.01459     0.00000
    0.00000     0.00000     0.00000     0.00000     0.00000     0.01459
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Voigt bulk  modulus, B_V =    89.33  GPa
    Voigt shear modulus, G_V =    68.61  GPa

    Reuss bulk  modulus, B_R =    89.33  GPa
    Reuss shear modulus, G_R =    68.61  GPa

    Hill bulk  modulus,  B_H =    89.33  GPa
    Hill shear modulus,  G_H =    68.61  GPa

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Voigt Young modulus,  E_V =   163.87  GPa
    Voigt Poisson ratio, nu_V =     0.19

    Reuss Young modulus,  E_R =   163.87  GPa
    Reuss Poisson ratio, nu_R =     0.19

    Hill Young modulus,   E_H =   163.87  GPa
    Hill Poisson ratio,  nu_H =     0.19

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Elastic Anisotropy in polycrystalline, AVR =    0.000 %
    Eigenvalues of elastic constant (stiffness) matrix:

           268.0
           137.4
           137.4
            68.5
            68.5
            68.5

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pughs modulus ratio, k=G/B (Voigt):    0.768
    Pughs modulus ratio, k=G/B (Reuss):    0.768
    Pughs modulus ratio, k=G/B (Hill) :    0.768

    Note:
    The Pughs modulus ratio k equals shear modulus divided by bulk modulus.
    It is used to estimate ductility and brittleness of materials.
    A value below 0.5 suggests ductile behavior, while above 0.5 indicates brittleness.
    This ratio is useful in mechanical design and failure prediction.
    It provides insight into bonding nature and structural integrity.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Vickers hardnesses (Voigt) (Tian  model, 2012):   13.601 [GPa]
    Vickers hardnesses (Voigt) (Chen  model, 2011):   14.425 [GPa]
    Vickers hardnesses (Voigt) (Teter model, 1998):   10.359 [GPa]
    Vickers hardnesses (Reuss) (Tian  model, 2012):   13.601 [GPa]
    Vickers hardnesses (Reuss) (Chen  model, 2011):   14.425 [GPa]
    Vickers hardnesses (Reuss) (Teter model, 1998):   10.359 [GPa]
    Vickers hardnesses (Hill)  (Tian  model, 2012):   13.601 [GPa]
    Vickers hardnesses (Hill)  (Chen  model, 2011):   14.425 [GPa]
    Vickers hardnesses (Hill)  (Teter model, 1998):   10.359 [GPa]

    Note (VASP case)
    Hill = Voigt-Reuss-Hill (VRH) Approximation (averages)
    For HvChen (= Vickers hardnesses (Chen  model, 2011))
    covalent and ionic crystals: Root Mean Square Error (RMSE) = 4.4
    covalent and ionic crystals: Mean Absolute Error (MAE)     = 2.1
    bulk metallic glasses:       Root Mean Square Error (RMSE) = 0.9
    bulk metallic glasses:       Mean Absolute Error (MAE)     = 0.8

    Note (Vickers Hardness and Structural Implications)
    Vickers hardness Hv estimates resistance to plastic deformation.
    It correlates with yield strength and fracture toughness.
    Useful for predicting wear resistance and indentation behavior.
    High Hv indicates strong atomic bonding and low ductility.
    Hv helps assess crack initiation and propagation under stress.
    Important for evaluating surface stability and coating durability.
    Guides material selection for cutting tools and structural components.
    Also used to estimate void formation resistance and microstructural integrity.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Universal Anisotropy Index:    0.000

    Note: Universal Elastic Anisotropy Index and Microstructural Effects
    Au measures overall elastic anisotropy in polycrystalline materials.
    Influences GP zone formation by controlling stress field directionality.
    Helps assess delamination risk and interfacial energy variation.
    Guides stacking orientation in composite design for stress optimization.
    Useful for predicting void formation in anisotropic stress environments.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Zener Anisotropy Ratio:    0.998

    Note: Zener Ratio and Structural Behavior
    Zener ratio quantifies elastic anisotropy in cubic crystals.
    It helps predict crack propagation direction in anisotropic stress fields.
    Guides selection of slip systems for dislocation motion.
    Important for epitaxial growth where lattice matching affects stability.
    Also used to evaluate directional stress concentration at interfaces.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Gruneisen constant (Voigt) (Slater approximation):    1.121
    Gruneisen constant (Reuss) (Slater approximation):    1.121
    Gruneisen constant (Hill)  (Slater approximation):    1.121

    Note:
    The Gruneisen parameter gamma describes lattice anharmonicity.
    In the Slater approximation, gamma is estimated from bulk modulus B and shear modulus G.
    This method does not require phonon data and enables fast thermal property estimation.
    Gamma affects the following properties:
      - Thermal expansion coefficient: alpha = gamma * Cv / (B * V)
      - Heat capacity temperature dependence (Debye model correction)
      - Lattice thermal conductivity (via phonon-phonon scattering, e.g., Slack model)
      - Phonon lifetime and scattering rate (related to anharmonicity)
      - Pressure dependence of vibrational properties and equation of state
      - Melting point under pressure (Lindemann theory)
      - Temperature dependence of elastic constants
    Even as an approximation, gamma links mechanical stiffness to thermal behavior.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vl = sound velocity (longitudinal wave) [m/s], vt = sound velocity (transverse waves) [m/s]
    vl/vt ratio:    1.623

    Note (Hill Average and Physical Property Relations)
    Hill average is the mean of Voigt and Reuss bounds for bulk and shear moduli.
    It provides stable input for property models derived from elastic constants.
    Used in Gruneisen parameter estimation via Slater approximation:
        gamma is calculated from bulk modulus B and shear modulus G.
    Used in hardness models:
        Hv is empirically modeled as a function of G and B in Tian, Chen, and Teter formulas.
    Used in ductility prediction: Pugh ratio k = G/B.
    Used in sound velocity estimation: vl = sqrt((B + 4/3 G)/rho), vt = sqrt(G/rho).
    vl/vt ratio relates to Poisson ratio and elastic anisotropy.
    These relations support thermal conductivity, EOS, and mechanical stability analysis.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    ... Have a G00D Day, Week, Month, Year, and Century (if you are lucky) ...
               Bye-Bye! Tschuess! Ciao! Poka! Zia Jian! KhodaHafez!

see ElaStic_2nd.out file in Energy-vs-Strain directory.
```

## Results: Stress and 2nd order: Diamond Si (Time: about 2 [min])
```
ElaStic_2nd.out found in Stress-vs-Strain. Displaying contents:
    The output of ElaStic code
    Today is Sat Aug 16 00:51:41 2025

    Symmetry of the second-order elastic constant matrix in Voigt notation.
    for, space group-number between 207 and 230, Cubic I structure.

               C11     C12     C12      0       0       0
               C12     C11     C12      0       0       0
               C12     C12     C11      0       0       0
                0       0       0      C44      0       0
                0       0       0       0      C44      0
                0       0       0       0       0      C44

    Elastic constant (stiffness) matrix in GPa:


      165.2        43.8        43.8        -0.0        -0.0        -0.0
       43.8       165.2        43.8        -0.0        -0.0        -0.0
       43.8        43.8       165.2        -0.0        -0.0        -0.0
       -0.0        -0.0        -0.0        67.8        -0.0        -0.0
       -0.0        -0.0        -0.0        -0.0        67.8        -0.0
       -0.0        -0.0        -0.0        -0.0        -0.0        67.8

    Elastic compliance matrix in 1/GPa:


    0.00681    -0.00143    -0.00143     0.00000     0.00000     0.00000
   -0.00143     0.00681    -0.00143     0.00000     0.00000     0.00000
   -0.00143    -0.00143     0.00681     0.00000     0.00000     0.00000
    0.00000     0.00000     0.00000     0.01474     0.00000     0.00000
    0.00000     0.00000     0.00000     0.00000     0.01474     0.00000
    0.00000     0.00000     0.00000     0.00000     0.00000     0.01474

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Voigt bulk  modulus, B_V =    84.28  GPa
    Voigt shear modulus, G_V =    64.97  GPa

    Reuss bulk  modulus, B_R =    84.28  GPa
    Reuss shear modulus, G_R =    64.78  GPa

    Hill bulk  modulus,  B_H =    84.28  GPa
    Hill shear modulus,  G_H =    64.87  GPa

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Voigt Young modulus,  E_V =   155.06  GPa
    Voigt Poisson ratio, nu_V =     0.19

    Reuss Young modulus,  E_R =   154.69  GPa
    Reuss Poisson ratio, nu_R =     0.19

    Hill Young modulus,   E_H =   154.88  GPa
    Hill Poisson ratio,  nu_H =     0.19

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Elastic Anisotropy in polycrystalline, AVR =    0.150 %

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Eigenvalues of elastic constant (stiffness) matrix:

           121.3
           252.8
           121.3
            67.8
            67.8
            67.8

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pughs modulus ratio, k=G/B (Voigt):    0.771
    Pughs modulus ratio, k=G/B (Reuss):    0.769
    Pughs modulus ratio, k=G/B (Hill) :    0.770

    Note:
    The Pughs modulus ratio k equals shear modulus divided by bulk modulus.
    It is used to estimate ductility and brittleness of materials.
    A value below 0.5 suggests ductile behavior, while above 0.5 indicates brittleness.
    This ratio is useful in mechanical design and failure prediction.
    It provides insight into bonding nature and structural integrity.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Vickers hardnesses (Voigt) (Tian  model, 2012):   13.143 [GPa]
    Vickers hardnesses (Voigt) (Chen  model, 2011):   13.953 [GPa]
    Vickers hardnesses (Voigt) (Teter model, 1998):    9.810 [GPa]
    Vickers hardnesses (Reuss) (Tian  model, 2012):   13.071 [GPa]
    Vickers hardnesses (Reuss) (Chen  model, 2011):   13.864 [GPa]
    Vickers hardnesses (Reuss) (Teter model, 1998):    9.781 [GPa]
    Vickers hardnesses (Hill)  (Tian  model, 2012):   13.107 [GPa]
    Vickers hardnesses (Hill)  (Chen  model, 2011):   13.909 [GPa]
    Vickers hardnesses (Hill)  (Teter model, 1998):    9.796 [GPa]

    Note (VASP case)
    Hill = Voigt-Reuss-Hill (VRH) Approximation (averages)
    For HvChen (= Vickers hardnesses (Chen  model, 2011))
    covalent and ionic crystals: Root Mean Square Error (RMSE) = 4.4
    covalent and ionic crystals: Mean Absolute Error (MAE)     = 2.1
    bulk metallic glasses:       Root Mean Square Error (RMSE) = 0.9
    bulk metallic glasses:       Mean Absolute Error (MAE)     = 0.8

    Note (Vickers Hardness and Structural Implications)
    Vickers hardness Hv estimates resistance to plastic deformation.
    It correlates with yield strength and fracture toughness.
    Useful for predicting wear resistance and indentation behavior.
    High Hv indicates strong atomic bonding and low ductility.
    Hv helps assess crack initiation and propagation under stress.
    Important for evaluating surface stability and coating durability.
    Guides material selection for cutting tools and structural components.
    Also used to estimate void formation resistance and microstructural integrity.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Universal Anisotropy Index:    0.015

    Note: Universal Elastic Anisotropy Index and Microstructural Effects
    Au measures overall elastic anisotropy in polycrystalline materials.
    Influences GP zone formation by controlling stress field directionality.
    Helps assess delamination risk and interfacial energy variation.
    Guides stacking orientation in composite design for stress optimization.
    Useful for predicting void formation in anisotropic stress environments.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Zener Anisotropy Ratio:    1.118

    Note: Zener Ratio and Structural Behavior
    Zener ratio quantifies elastic anisotropy in cubic crystals.
    It helps predict crack propagation direction in anisotropic stress fields.
    Guides selection of slip systems for dislocation motion.
    Important for epitaxial growth where lattice matching affects stability.
    Also used to evaluate directional stress concentration at interfaces.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Gruneisen constant (Voigt) (Slater approximation):    1.120
    Gruneisen constant (Reuss) (Slater approximation):    1.120
    Gruneisen constant (Hill)  (Slater approximation):    1.120

    Note:
    The Gruneisen parameter gamma describes lattice anharmonicity.
    In the Slater approximation, gamma is estimated from bulk modulus B and shear modulus G.
    This method does not require phonon data and enables fast thermal property estimation.
    Gamma affects the following properties:
      - Thermal expansion coefficient: alpha = gamma * Cv / (B * V)
      - Heat capacity temperature dependence (Debye model correction)
      - Lattice thermal conductivity (via phonon-phonon scattering, e.g., Slack model)
      - Phonon lifetime and scattering rate (related to anharmonicity)
      - Pressure dependence of vibrational properties and equation of state
      - Melting point under pressure (Lindemann theory)
      - Temperature dependence of elastic constants
    Even as an approximation, gamma links mechanical stiffness to thermal behavior.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vl = sound velocity (longitudinal wave) [m/s], vt = sound velocity (transverse waves) [m/s]
    vl/vt ratio:    1.622

    Note (Hill Average and Physical Property Relations)
    Hill average is the mean of Voigt and Reuss bounds for bulk and shear moduli.
    It provides stable input for property models derived from elastic constants.
    Used in Gruneisen parameter estimation via Slater approximation:
        gamma is calculated from bulk modulus B and shear modulus G.
    Used in hardness models:
        Hv is empirically modeled as a function of G and B in Tian, Chen, and Teter formulas.
    Used in ductility prediction: Pugh ratio k = G/B.
    Used in sound velocity estimation: vl = sqrt((B + 4/3 G)/rho), vt = sqrt(G/rho).
    vl/vt ratio relates to Poisson ratio and elastic anisotropy.
    These relations support thermal conductivity, EOS, and mechanical stability analysis.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    ... Have a G00D Day, Week, Month, Year, and Century (if you are lucky) ...
               Bye-Bye! Tschuess! Ciao! Poka! Zia Jian! KhodaHafez!
```

## Building data in a Python environment
1. sudo apt -y install python3-pip
2. pip freeze > requirements.txt
