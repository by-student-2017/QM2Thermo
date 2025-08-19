#!/bin/bash
#
EXECUTABLE=$EXCITINGROOT/bin/exciting_purempi

# Set number of threads and CPUs
export OMP_NUM_THREADS=1
NCPUs=$(($(nproc) / 2))

label=`ls -d Dst??`
for Dstn in $label ; do
    cd $Dstn
    Dstn_num_list=`ls -d ${Dstn}_??`
    for Dstn_num in $Dstn_num_list ; do
        cd $Dstn_num/
        cp -f $Dstn_num.xml input.xml
        echo
        echo '        +--------------------------------------+'
        echo '        | SCF calculation of "'${Dstn_num}'" starts |'
        echo '        +--------------------------------------+'
        time mpirun -np ${NCPUs} $EXECUTABLE | tee output.screen
        date
        cd ../
    done
    cd ../
done
