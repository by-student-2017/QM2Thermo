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

## Result: Diamond C
- (Ref: C11 = 1076, C12 = 125, C44 = 577, B0 = 452 [GPa])
```
ElaStic_2nd.out found in Energy-vs-Strain. Displaying contents:
    The output of ElaStic code
    Today is Mon Aug 18 15:22:18 2025

    Symmetry of the second-order elastic constant matrix in Voigt notation.
    for, space group-number between 207 and 230, Cubic I structure.

               C11     C12     C12      0       0       0
               C12     C11     C12      0       0       0
               C12     C12     C11      0       0       0
                0       0       0      C44      0       0
                0       0       0       0      C44      0
                0       0       0       0       0      C44

    Elastic constant (stiffness) matrix in GPa:

      1079.1       134.2       134.2         0.0         0.0         0.0
       134.2      1079.1       134.2         0.0         0.0         0.0
       134.2       134.2      1079.1         0.0         0.0         0.0
         0.0         0.0         0.0       571.7         0.0         0.0
         0.0         0.0         0.0         0.0       571.7         0.0
         0.0         0.0         0.0         0.0         0.0       571.7


    Elastic compliance matrix in 1/GPa:

     0.00095    -0.00011    -0.00011     0.00000     0.00000     0.00000
    -0.00011     0.00095    -0.00011     0.00000     0.00000     0.00000
    -0.00011    -0.00011     0.00095     0.00000     0.00000     0.00000
     0.00000     0.00000     0.00000     0.00175     0.00000     0.00000
     0.00000     0.00000     0.00000     0.00000     0.00175     0.00000
     0.00000     0.00000     0.00000     0.00000     0.00000     0.00175

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Voigt bulk  modulus, B_V =   449.17  GPa
    Voigt shear modulus, G_V =   531.98  GPa

    Reuss bulk  modulus, B_R =   449.17  GPa
    Reuss shear modulus, G_R =   527.36  GPa

    Hill bulk  modulus,  B_H =   449.17  GPa
    Hill shear modulus,  G_H =   529.67  GPa

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Voigt Young modulus,  E_V =  1144.22  GPa
    Voigt Poisson ratio, nu_V =     0.08

    Reuss Young modulus,  E_R =  1137.08  GPa
    Reuss Poisson ratio, nu_R =     0.08

    Hill Young modulus,   E_H =  1140.65  GPa
    Hill Poisson ratio,  nu_H =     0.08

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Elastic Anisotropy in polycrystalline, AVR =    0.436 %

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Eigenvalues of elastic constant (stiffness) matrix:

           944.8
          1347.5
           944.8
           571.7
           571.7
           571.7

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pughs modulus ratio, k=G/B (Voigt):    1.184
    Pughs modulus ratio, k=G/B (Reuss):    1.174
    Pughs modulus ratio, k=G/B (Hill) :    1.179

    Note:
    The Pughs modulus ratio k equals shear modulus divided by bulk modulus.
    It is used to estimate ductility and brittleness of materials.
    A value below 0.5 suggests ductile behavior, while above 0.5 indicates brittleness.
    This ratio is useful in mechanical design and failure prediction.
    It provides insight into bonding nature and structural integrity.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Vickers hardnesses (Voigt) (Tian  model, 2012):   94.899 [GPa]
    Vickers hardnesses (Voigt) (Chen  model, 2011):   92.862 [GPa]
    Vickers hardnesses (Voigt) (Teter model, 1998):   80.328 [GPa]
    Vickers hardnesses (Reuss) (Tian  model, 2012):   93.385 [GPa]
    Vickers hardnesses (Reuss) (Chen  model, 2011):   91.406 [GPa]
    Vickers hardnesses (Reuss) (Teter model, 1998):   79.631 [GPa]
    Vickers hardnesses (Hill)  (Tian  model, 2012):   94.141 [GPa]
    Vickers hardnesses (Hill)  (Chen  model, 2011):   92.133 [GPa]
    Vickers hardnesses (Hill)  (Teter model, 1998):   79.980 [GPa]

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

    Universal Anisotropy Index:    0.044

    Note: Universal Elastic Anisotropy Index and Microstructural Effects
    Au measures overall elastic anisotropy in polycrystalline materials.
    Influences GP zone formation by controlling stress field directionality.
    Helps assess delamination risk and interfacial energy variation.
    Guides stacking orientation in composite design for stress optimization.
    Useful for predicting void formation in anisotropic stress environments.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Zener Anisotropy Ratio:    1.210

    Note: Zener Ratio and Structural Behavior
    Zener ratio quantifies elastic anisotropy in cubic crystals.
    It helps predict crack propagation direction in anisotropic stress fields.
    Guides selection of slip systems for dislocation motion.
    Important for epitaxial growth where lattice matching affects stability.
    Also used to evaluate directional stress concentration at interfaces.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Gruneisen constant (Voigt) (Slater approximation):    1.041
    Gruneisen constant (Reuss) (Slater approximation):    1.042
    Gruneisen constant (Hill)  (Slater approximation):    1.042

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
    vl/vt ratio:    1.477

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



