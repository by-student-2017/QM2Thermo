#!/bin/bash

filename=wien
cp ./$filename/*.dos1    wien.dos1
cp ./$filename/*.energy  wien.energy
cp ./$filename/*.kgen    wien.kgen
cp ./$filename/*.klist   wien.klist
cp ./$filename/*.struct  wien.struct

if [ ! -s "cfA1.dat" ]; then
    echo "cfA1.dat is not exists. Running asupersi.exe."
    echo -e "\n ----- run generate_stencil.exe ----- \n"
    ./generate_stencil.exe
    echo -e "\n\n ----- run group_velocity.exe ------- \n"
    ./group_velocity.exe
    echo -e "\n\n ----- run chemical_potential.exe --- \n"
    ./chemical_potential.exe
    gnuplot plot_cp.gpl
    echo -e "\n\n ----- run Seebeck_analysis.exe ----- "
    ./seebeck_analysis.exe
    echo -e "\n\n ----- gnuplot ---------------------- \n"
    gnuplot plot_Seebeck.gpl
    gnuplot plot_ABGV2D.gpl
elif [ ! -s "AKK.DATA" ]; then
    echo "AKK.DATA is not exists. skip generate_stencil.exe"
    echo -e "\n ----- run group_velocity.exe ------- \n"
    ./group_velocity.exe
    echo -e "\n\n ----- run chemical_potential.exe --- \n"
    ./chemical_potential.exe
    gnuplot plot_cp.gpl
    echo -e "\n\n ----- run Seebeck_analysis.exe ----- "
    ./seebeck_analysis.exe
    echo -e "\n\n ----- gnuplot ---------------------- \n"
    gnuplot plot_Seebeck.gpl
    gnuplot plot_ABGV2D.gpl
elif [ ! -s "apot.data" ]; then
    echo "AKK.DATA exists. skip generate_stencil.exe"
    echo "apot.data not found. skip chemical_potential.exe"
    echo -e "\n ----- run chemical_potential.exe --- \n"
    ./chemical_potential.exe
    gnuplot plot_cp.gpl
    echo -e "\n\n ----- run Seebeck_analysis.exe ----- "
    ./seebeck_analysis.exe
    echo -e "\n\n ----- gnuplot ---------------------- \n"
    gnuplot plot_Seebeck.gpl
    gnuplot plot_ABGV2D.gpl
else
    echo "AKK.DATA exists. skip generate_stencil.exe"
    echo "apot.data exists. skip chemical_potential.exe"
    echo -e "\n ----- run Seebeck_analysis.exe ----- "
    ./seebeck_analysis.exe
    echo -e "\n ----- gnuplot ---------------------- \n"
    gnuplot plot_Seebeck.gpl
    gnuplot plot_ABGV2D.gpl
fi

echo "----------------------------------------------------------------------------------"
echo '## set terminal win font "Arial,12"                                             ##'
echo "## warning: unknown or ambiguous terminal type; type 'set terminal' for a list. ##"
echo ""
echo "The above 'warning' is for convenience on Windows, so don't worry about it."
echo "----------------------------------------------------------------------------------"
