# QE v7.2 + thermo_pw v1.8.1

## QE v7.2 Installation
1. sudo apt update
2. sudo apt -y install gfortran g++ build-essential make libopenblas-dev libopenmpi-dev libfftw3-dev
3. wget https://github.com/QEF/q-e/archive/refs/tags/qe-7.2.tar.gz
4. tar xvf qe-7.2.tar.gz
5. cd q-e-qe-7.2
6. ./configure
7. make pwall
8. sudo make install

## thermo_pw 1.8.1 Installation for QE 7.2
1. cd q-e-qe-7.2
2. wget https://github.com/dalcorso/thermo_pw/archive/refs/tags/1.8.1.tar.gz
3. tar -xvzf 1.8.1.tar.gz
4. mv thermo_pw-1.8.1 thermo_pw
5. cd thermo_pw
6. make join_qe
7. cd ..
8. ./configure
9. make thermo_pw
10. sudo make install

## Elastic constants calculation
1. bash ./run_elastic_thermo_pw.sh
2. grep -A 50 "Elastic" log/case.elastic.out | tail -70

## Restart
- You need to delete directories such as restart.

## If there is not enough memory
- Reduce the number of parallel processes as follows (in this example it is 1, but if you have enough memory you can use 2 or 4):
```
NCPUs=1
```

## Show results
- grep -A 50 "Elastic" case.elastic.out | tail -70
```
     Elastic constants C_ij (kbar)
    i j=        1           2           3           4           5           6
    1  1543.12415   567.35371   567.35371     0.00000     0.00000     0.00000
    2   567.35371  1543.12415   567.35371     0.00000     0.00000     0.00000
    3   567.35371   567.35371  1543.12415     0.00000     0.00000     0.00000
    4     0.00000     0.00000     0.00000   766.93681     0.00000     0.00000
    5     0.00000     0.00000     0.00000     0.00000   766.93681     0.00000
    6     0.00000     0.00000     0.00000     0.00000     0.00000   766.93681

     1 bar = 10^5 Pa; 10 kbar = 1 GPa; 1 atm = 1.01325 bar; 1 Pa = 1 N/m^2
     1 Pa = 10 dyn/cm^2; 1 Mbar = 10^11 Pa
     1 torr = 1 mm Hg = 1/760 bar = 7.5006 x 10^-3 Pa


                    ----------------------------------------


     Elastic compliances  S_ij (1/Mbar)
    i j=        1           2           3           4           5           6
    1     0.80770    -0.21713    -0.21713    -0.00000     0.00000    -0.00000
    2    -0.21713     0.80770    -0.21713    -0.00000     0.00000    -0.00000
    3    -0.21713    -0.21713     0.80770    -0.00000     0.00000    -0.00000
    4     0.00000     0.00000     0.00000     1.30389     0.00000    -0.00000
    5     0.00000     0.00000     0.00000     0.00000     1.30389    -0.00000
    6     0.00000     0.00000     0.00000     0.00000     0.00000     1.30389

     1/Mbar = 1/10^{11} Pa; 1 Pa = 1 N/m^2

                    ----------------------------------------


     Voigt approximation:
     Bulk modulus  B =    892.61052 kbar
     Young modulus E =   1579.43164 kbar
     Shear modulus G =    655.31617 kbar
     Poisson Ratio n =      0.20509

     Reuss approximation:
     Bulk modulus  B =    892.61052 kbar
     Young modulus E =   1518.49951 kbar
     Shear modulus G =    624.14258 kbar
     Poisson Ratio n =      0.21647

     Voigt-Reuss-Hill average of the two approximations:
     Bulk modulus  B =    892.61052 kbar
     Young modulus E =   1548.96558 kbar
     Shear modulus G =    639.72938 kbar
     Poisson Ratio n =      0.21064

     Voigt-Reuss-Hill average; sound velocities:

     Compressional V_P =     8747.682 m/s
     Bulk          V_B =     6255.380 m/s
     Shear         V_G =     5295.670 m/s

     The approximate Debye temperature is      637.299 K

                    ----------------------------------------

     Average Debye sound velocity =     5827.179 m/s

     Debye temperature =      634.507 K
```

# version information
| thermo_pw | QE    |
| 1.8.1     | 7.2   |
| 2.0.3     | 7.4.1 |

## Test
- OS: Ubuntu 22.04 LTS (WSL2, Windows 11)
- GPU: 12th Gen Intel(R) Core(TM) i7-12700
- Memory: 32 GB


# (uninstall thermo_pw)
1. cd q-e-qe-7.2
2. cd thermo_pw
3. make leave_qe
