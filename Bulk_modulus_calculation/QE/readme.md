# QE v6.8 Installation (Ubuntu 22.04 LTS)
1. sudo apt update
2. sudo apt -y install gfortran g++ build-essential make libopenblas-dev libopenmpi-dev libfftw3-dev
3. wget https://github.com/QEF/q-e/archive/refs/tags/qe-6.8.tar.gz
4. tar xvf qe-6.8.tar.gz
5. cd q-e-qe-6.8
6. ./configure
7. make pwall
8. sudo make install

# Shear modulus calculation
1. bash ./run_bulk_qe.sh
2. awk -f compute_bulk_modulus_from_stress_qe.awk shear_results.txt
