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
1. bash ./run_elastic_qe.sh
2. awk -f compute_elastic_constants_from_stress_qe.awk elastic_results.txt
3. python3 compliance_python3.py