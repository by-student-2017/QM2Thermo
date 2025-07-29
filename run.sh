#!/bin/bash

filename=wien
cp ./$filename/*.dos1    wien.dos1
cp ./$filename/*.energy  wien.energy
cp ./$filename/*.kgen    wien.kgen
cp ./$filename/*.klist   wien.klist
cp ./$filename/*.struct  wien.struct

if [ ! -s "cfA1.dat" ]; then
    echo "cfA1.dat is not exists. Running asupersi.exe."
    ./generate_stencil.exe
    ./group_velocity.exe
    ./chemical_potential.exe
    ./plot.gpl
    ./seebeck_analysis.exe
    ./plot_Seebeck.gpl
elif [ ! -s "AKK.DATA" ]; then
    echo "AKK.DATA is not exists. Running asupersi.exe."
    ./group_velocity.exe
    ./chemical_potential.exe
    ./plot.gpl
    ./seebeck_analysis.exe
    ./plot_Seebeck.gpl
elif [ ! -s "apot.data" ]; then
    echo "AKK.DATA exists. Running asupersi.exe."
    echo "apot.data not found. Running aaacp.exe instead."
    ./chemical_potential.exe
    ./plot_cp.gpl
    ./seebeck_analysis.exe
    ./plot_Seebeck.gpl
else
    echo "AKK.DATA exists. Running asupersi.exe."
    echo "apot.data exists. Running aaass1.exe instead."
    ./seebeck_analysis.exe
    ./plot_Seebeck.gpl
fi
