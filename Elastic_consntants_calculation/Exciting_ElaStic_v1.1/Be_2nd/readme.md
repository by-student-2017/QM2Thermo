# Usage

## Easy method
1. bash run_exciting_energy_2nd.sh

## Step by step
1. ElaStic_Setup
2. (input your calculation conditions)
3. bash ElaStic@exciting-submit.sh
4. ElaStic_exciting_analyze.py
   - or ElaStic_Analyze
5. ElaStic_Result
6. cat ElaStic_2nd.out

## Troubleshooting
- If it doesn't work, give it permission to run:
```
cd $HOME/ElaStic
chmod +x *
```
- Problems may occur if you are using Elastic2020 etc. at the same time. It is recommended to specify the absolute path and run it. An example is shown below.
```
python3 $HOME/ElaStic/ElaStic_Setup
# (input your calculation conditions)
bash $HOME/ElaStic/ElaStic@exciting-submit.sh
python3 $HOME/ElaStic/ElaStic_exciting_analyze.py
# or python3 $HOME/ElaStic/ElaStic_Analyze
python3 $HOME/ElaStic/ElaStic_Result
cat ElaStic_2nd.out
```

## Results: HCP Be
```
ElaStic_2nd.out found in Energy-vs-Strain. Displaying contents:
    The output of ElaStic code
    Today is Mon Aug 18 15:15:15 2025

    Symmetry of the second-order elastic constant matrix in Voigt notation.
    for, space group-number between 177 and 194, Hexagonal I structure.

               C11     C12     C13      0       0       0
               C12     C11     C13      0       0       0
               C13     C13     C33      0       0       0
                0       0       0      C44      0       0
                0       0       0       0      C44      0
                0       0       0       0       0   (C11-C12)/2

    Elastic constant (stiffness) matrix in GPa:

       357.2       -20.4         5.1         0.0         0.0         0.0
       -20.4       357.2         5.1         0.0         0.0         0.0
         5.1         5.1       416.3         0.0         0.0         0.0
         0.0         0.0         0.0       178.8         0.0         0.0
         0.0         0.0         0.0         0.0       178.8         0.0
         0.0         0.0         0.0         0.0         0.0       188.8


    Elastic compliance matrix in 1/GPa:

     0.00281     0.00016    -0.00004     0.00000     0.00000     0.00000
     0.00016     0.00281    -0.00004     0.00000     0.00000     0.00000
    -0.00004    -0.00004     0.00240     0.00000     0.00000     0.00000
     0.00000     0.00000     0.00000     0.00559     0.00000     0.00000
     0.00000     0.00000     0.00000     0.00000     0.00559     0.00000
     0.00000     0.00000     0.00000     0.00000     0.00000     0.00530

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Voigt bulk  modulus, B_V =   123.36  GPa
    Voigt shear modulus, G_V =   185.33  GPa

    Reuss bulk  modulus, B_R =   121.97  GPa
    Reuss shear modulus, G_R =   184.76  GPa

    Hill bulk  modulus,  B_H =   122.67  GPa
    Hill shear modulus,  G_H =   185.05  GPa

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Voigt Young modulus,  E_V =   370.46  GPa
    Voigt Poisson ratio, nu_V =    -0.00

    Reuss Young modulus,  E_R =   368.32  GPa
    Reuss Poisson ratio, nu_R =    -0.00

    Hill Young modulus,   E_H =   369.39  GPa
    Hill Poisson ratio,  nu_H =    -0.00

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Elastic Anisotropy in polycrystalline, AVR =    0.153 %

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Eigenvalues of elastic constant (stiffness) matrix:

           336.2
           377.6
           416.9
           178.8
           178.8
           188.8

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pughs modulus ratio, k=G/B (Voigt):    1.502
    Pughs modulus ratio, k=G/B (Reuss):    1.515
    Pughs modulus ratio, k=G/B (Hill) :    1.509

    Note:
    The Pughs modulus ratio k equals shear modulus divided by bulk modulus.
    It is used to estimate ductility and brittleness of materials.
    A value below 0.5 suggests ductile behavior, while above 0.5 indicates brittleness.
    This ratio is useful in mechanical design and failure prediction.
    It provides insight into bonding nature and structural integrity.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Vickers hardnesses (Voigt) (Tian  model, 2012):   58.951 [GPa]
    Vickers hardnesses (Voigt) (Chen  model, 2011):   65.330 [GPa]
    Vickers hardnesses (Voigt) (Teter model, 1998):   27.985 [GPa]
    Vickers hardnesses (Reuss) (Tian  model, 2012):   59.375 [GPa]
    Vickers hardnesses (Reuss) (Chen  model, 2011):   65.865 [GPa]
    Vickers hardnesses (Reuss) (Teter model, 1998):   27.899 [GPa]
    Vickers hardnesses (Hill)  (Tian  model, 2012):   59.162 [GPa]
    Vickers hardnesses (Hill)  (Chen  model, 2011):   65.596 [GPa]
    Vickers hardnesses (Hill)  (Teter model, 1998):   27.942 [GPa]

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

    Universal Anisotropy Index:    0.027

    Note: Universal Elastic Anisotropy Index and Microstructural Effects
    Au measures overall elastic anisotropy in polycrystalline materials.
    Influences GP zone formation by controlling stress field directionality.
    Helps assess delamination risk and interfacial energy variation.
    Guides stacking orientation in composite design for stress optimization.
    Useful for predicting void formation in anisotropic stress environments.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Zener Anisotropy Ratio:    0.947

    Note: Zener Ratio and Structural Behavior
    Zener ratio quantifies elastic anisotropy in cubic crystals.
    It helps predict crack propagation direction in anisotropic stress fields.
    Guides selection of slip systems for dislocation motion.
    Important for epitaxial growth where lattice matching affects stability.
    Also used to evaluate directional stress concentration at interfaces.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Gruneisen constant (Voigt) (Slater approximation):    1.000
    Gruneisen constant (Reuss) (Slater approximation):    0.998
    Gruneisen constant (Hill)  (Slater approximation):    0.999

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
    vl/vt ratio:    1.413

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


