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

## Calculation method (stress base. 2nd order)
1. run_elastic_stress_2nd.sh

## Calculation method (energy base. 2nd order)
1. run_elastic_energy_2nd.sh

## Calculation method (stress base. 3rd order) (Some confirmation required)
1. run_elastic_stress_3rd.sh

## Calculation method (energy base. 3rd order) (All need to be checked)
1. run_elastic_energy_3rd.sh

## Building data in a Python environment
1. sudo apt -y install python3-pip
2. pip freeze > requirements.txt
