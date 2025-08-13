#!/bin/bash

# Base input file
base_input="case.scf.in" # old version

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
input_file=log/case.elastic.in
output_file=log/case.elastic.out

cp ${base_input} ${input_file}

## For ibrav != 0 case
#  what='elastic_constants_geo',
#  frozen_ions=.FALSE.,
#  lmurn=.FALSE.,

## For ibrav = 0 case
#  what='scf_elastic_constants',
#  frozen_ions=.FALSE.,
#  lmurn=.FALSE.,
#  find_ibrav=.TRUE.,

## For using Murnahan
#  what = 'mur_lc_elastic_constants',
#  frozen_ions = .FALSE.,
#  lmurn = .TRUE.,
#  find_ibrav = .TRUE.,
#  ngeo = 7,
#  pmin = -5.0,
#  pmax = 5.0,

cat << EOF > thermo_control
&INPUT_THERMO
  what='scf_elastic_constants',
  frozen_ions=.FALSE.,
  lmurn=.FALSE.,
  find_ibrav=.TRUE.,
/
EOF

# thermo_pw execution (elastic constant tensor)
${PARA_IMAGE_PREFIX} ${BIN_DIR}/thermo_pw.x < ${input_file} | tee ${output_file}

grep -A 50 "Elastic" ${output_file} | tail -70
