# QE v6.8 Installation (Ubuntu 22.04 LTS)
1. sudo apt update
2. sudo apt -y install gfortran g++ build-essential make libopenblas-dev libopenmpi-dev libfftw3-dev
3. wget https://github.com/QEF/q-e/archive/refs/tags/qe-6.8.tar.gz
4. tar xvf qe-6.8.tar.gz
5. cd q-e-qe-6.8
6. ./configure
7. make pwall
8. sudo make install

# python
1. pip install numpy

# Shear modulus calculation
1. bash ./run_opt.sh
2. bash ./run_elastic_qe.sh
3. awk -f compute_elastic_constants_from_stress_qe.awk elastic_results.txt
4. python3 compliance_python3.py

# If there is not enough memory
- Reduce the number of parallel processes as follows (in this example it is 1, but if you have enough memory you can use 2 or 4):
```
NCPUs=1
```

# Note
- For conventional cells, values close to those of Materials Projects were obtained. However, for primitive cells, large values were obtained for C14, C56, etc.
- It is better to use relax to optimize the atomic positions. For structures other than the reference structure, calculations will be performed with the settings unchanged, so you can specify relax in case.scf.in or case.opt.in.
- Si Data (elastic.txt) in this code. (Conventional cell calculation. see case.opt.in) (calculation time: 1 h)
```
=========================================
Components of the Elastic Constant Tensor
=========================================
Elastic Constant C11all =    153.93802725 [GPa]
Elastic Constant C22all =    153.98215875 [GPa]
Elastic Constant C33all =    153.96744825 [GPa]
Elastic Constant C12all =     56.77149712 [GPa]
Elastic Constant C13all =     56.76414187 [GPa]
Elastic Constant C23all =     56.79724050 [GPa]
Elastic Constant C44all =     74.80657013 [GPa]
Elastic Constant C55all =     74.81024775 [GPa]
Elastic Constant C66all =     74.79185963 [GPa]
Elastic Constant C14all =     -0.00183881 [GPa]
Elastic Constant C15all =      0.00183881 [GPa]
Elastic Constant C16all =      0.00000000 [GPa]
Elastic Constant C24all =      0.00000000 [GPa]
Elastic Constant C25all =      0.00551644 [GPa]
Elastic Constant C26all =      0.00367762 [GPa]
Elastic Constant C34all =      0.00367762 [GPa]
Elastic Constant C35all =      0.00183881 [GPa]
Elastic Constant C36all =     -0.00367762 [GPa]
Elastic Constant C45all =     -0.00551644 [GPa]
Elastic Constant C46all =     -0.00183881 [GPa]
Elastic Constant C56all =     -0.00367762 [GPa]
=========================================
Average properties for a cubic crystal
=========================================
Bulk Modulus    =     89.17259925 [GPa]
Shear Modulus 1 =     74.80289250 [GPa]
Shear Modulus 2 =     48.59245913 [GPa]
Poisson's ratio =      0.26942005 [GPa]
```

## Test
- OS: Ubuntu 22.04 LTS (WSL2, Windows 11)
- GPU: 12th Gen Intel(R) Core(TM) i7-12700
- Memory: 32 GB
