# Exciting neon-21 + ElaStic v1.1

## Exciting neon-21 Installation
1. sudo apt update
2. sudo apt -y install gfortran g++ build-essential make libopenblas-dev libopenmpi-dev libfftw3-dev
3. sudo apt -y install gnuplot ghostscript
4. sudo apt -y install python3-all-dev graphviz xsltproc
5. sudo apt install hdf5-helpers libhdf5-dev libhdf5-openmpi-dev fftw3-dev
6. pip3 install lxml
7. cd $HOME
8. mkdir exciting.neon-21
9. cd exciting.neon-21
10. wget https://exciting-code.org/uploads/exciting/exciting.neon-21.tar.gz
11. tar xvf exciting.neon-21.tar.gz
12. cp build/platforms/make.inc.gfortran10+.hdf5.fftw3 build/make.inc
13. make mpi
14. echo 'export EXCITINGROOT="$HOME/exciting.neon-21"' >> ~/.bashrc
15. echo 'export PATH="$EXCITINGROOT/bin:$PATH"' >> ~/.bashrc
16. echo 'export PATH="$EXCITINGROOT/tools/tutorial_scripts:$PATH"' >> ~/.bashrc

## Elastic Installation
1. sudo apt update
2. sudo apt -y install gfortran make build-essential grace
3. sudo apt -y install python3 python3-numpy python3-matplotlib
4. cd $HOME
5. git clone https://github.com/by-student-2017/ElaStic.git
6. cd ElaStic
7. tar xvf adon_v1_0.tar.gz
8. cd SpaceGroups
9. make
10. cd ..
11. cp SpaceGroups/sgroup ./sgroup
12. chmod +x *
13. echo 'export ElaSticROOT=$HOME/ElaStic' >> ~/.bashrc
14. echo 'export PATH="$ElaSticROOT:$PATH"' >> ~/.bashrc

## Tutorial
- The tutorial directory contains a readme.md file that explains how to calculate the code. Please refer to it.

## Test
- OS: Ubuntu 22.04 LTS (WLS2, Windows 11)
- Python: 3.10.12
- Numpy: 1.21.5
- Matplotlib: 3.5.1
- lxml: 4.8.0
- exciting: neon-21
- gfortran: gcc version 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04)

## References
- [1] R. Golesorkhtabar, P. Pavone, J. Spitaler, P. Puschnig, and C. Draxl, ElaStic: A tool for calculating second-order elastic constants from first principles, Comp. Phys. Commun. 184, 1861 (2013).: https://exciting-code.org/elastic/
- [2] excitingscripts: https://github.com/sr76/excitingscripts/
- [3] excitingscripts 1.0.0: https://pypi.org/project/excitingscripts/
- [4] Elastic1.0-wien2k: https://github.com/jazairdz/Elastic1.0-wien2k
- [5] WT16: https://rashid-phy.github.io/me/wien2k.html
- [6] https://www.youtube.com/@phy22: https://www.youtube.com/hashtag/physicsschool20
- [7] exciting-neon.21/tools/tutorial_scripts: http://exciting.wikidot.com/neon-how-to-calculate-stress-tensor
