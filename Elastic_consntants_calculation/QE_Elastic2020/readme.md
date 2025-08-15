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

## Calculation method (energy base. 2nd order)
1. run_elastic_energy_2nd.sh

## Calculation method (stress base. 2nd order)
1. run_elastic_stress_2nd.sh

## Calculation method (stress base. 3rd order) (Some confirmation required: need to check RII to N structure)
1. run_elastic_stress_3rd.sh

## Calculation method (energy base. 3rd order) (All need to be checked)
1. run_elastic_energy_3rd.sh

## Whether or not to display the diagram
- The figure is stored in a directory named Plot. If you want to see it immediately (*.sh), do the following:
```
# Choose whether to display the diagram: yes or no
show_plot_flag="yes"
```

## Test
- OS: Ubuntu 22.04 LTS (WLS2)
- Python: 3.10.12
- Numpy: 1.21.5
- Matplotlib: 3.5.1

## Results: Energy and 2nd order
```
ElaStic_2nd.out found in Energy-vs-Strain. Displaying contents:
    The output of ElaStic code
    Today is Fri Aug 15 13:29:00 2025

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

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Eigenvalues of elastic constant (stiffness) matrix:

           268.0
           137.4
           137.4
            68.5
            68.5
            68.5

    ... Have a G00D Day, Week, Month, Year, and Century (if you are lucky) ...
               Bye-Bye! Tschuess! Ciao! Poka! Zia Jian! KhodaHafez!
```

## Results: Stress and 2nd order
```
ElaStic_2nd.out found in Stress-vs-Strain. Displaying contents:
    The output of ElaStic code
    Today is Fri Aug 15 13:30:19 2025

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

    ... Have a G00D Day, Week, Month, Year, and Century (if you are lucky) ...
               Bye-Bye! Tschuess! Ciao! Poka! Zia Jian! KhodaHafez!
```

## Building data in a Python environment
1. sudo apt -y install python3-pip
2. pip freeze > requirements.txt
