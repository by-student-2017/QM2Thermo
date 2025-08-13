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
1. cd $HOME
2. git clone https://github.com/ponychen123/Elastic2020.git
3. cd Elastic2020
4. python3 -m venv Elastic2020
5. source Elastic2020/bin/activate
6. pip3 install numpy==1.23 matplotlib
7. tar xvf adon_v1_0.tar.gz
8. cd SpaceGroups
9. make
10. cd ..
11. cp SpaceGroups/sgroup ./sgroup
12. echo 'export PATH="$HOME/Elastic2020:$PATH"' >> ~/.bashrc

## Calculation method
1. bash run_elastic.sh
