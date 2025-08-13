#!/bin/bash

# Base input file
base_input="case_ibrav.scf.in" # old version

# Set number of threads and CPUs
export OMP_NUM_THREADS=1
NCPUs=$(($(nproc) / 2))

# Environment variables (modify as needed)
BIN_DIR=/usr/local/bin
TMP_DIR=./work
PARA_IMAGE_PREFIX="mpirun -np ${NCPUs}"

# make work directory
mkdir -p work
mkdir -p log

# input file
input_file="log/case.disp.in"
output_file="log/case.disp.out"

cp ${base_input} ${input_file}

#  find_ibrav=.TRUE.,
#  continue_zero_ibrav=.TRUE.,

cat > thermo_control << EOF
 &INPUT_THERMO
  what='mur_lc_disp',
  nq1_d=128, nq2_d=128, nq3_d=128,
 /
EOF

cat > ph_control << EOF
# Comment line: phonon calculation
&inputph
  tr2_ph = 1.0d-12 ,
  prefix='POSCAR',
  outdir='./work/g1/',
  fildyn = 'case.dyn.xml' ,
  ldisp = .TRUE. ,
  nq1=2, nq2=2, nq3=2,
/
EOF

# thermo_pw execution (dispersion)
${PARA_IMAGE_PREFIX} ${BIN_DIR}/thermo_pw.x  < ${input_file} | tee ${output_file}

ps2pdf output_disp.g1.ps disp.pdf
ps2pdf output_dos.g1.ps dos.pdf
ps2pdf output_therm.g1.ps therm.pdf