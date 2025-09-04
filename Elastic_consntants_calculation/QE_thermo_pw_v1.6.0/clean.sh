#!/bin/bash

echo "# delet old files"
echo ""

echo "delet output files from run_elastic_thermo_pw.sh"
rm -rf log work restart therm_files gnuplot_files elastic_constants thermo_control *.ps
echo ""

echo "delet output files from run_disp_thermo_pw.sh"
rm -rf log work restart therm_files gnuplot_files thermo_control *.ps
rm -rf phdisp_files dynamical_matrices ph_control *.pdf
echo ""

echo "delet output files from run_disp_mur_thermo_pw.sh"
rm -rf log work restart therm_files gnuplot_files thermo_control *.ps
rm -rf phdisp_files dynamical_matrices ph_control *.pdf
rm -rf energy_files
echo ""

echo "delet output files from run_Gruneisen_thermo_pw.sh"
rm -rf log work restart therm_files gnuplot_files elastic_constants thermo_control *.ps
rm -rf phdisp_files dynamical_matrices ph_control *.pdf
rm -rf energy_files
rm -rf anhar_files
