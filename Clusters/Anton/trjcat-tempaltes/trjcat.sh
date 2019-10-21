#!/bin/bash

nn=0
for ii in $(seq 0 10000 90000)
do
    nn=$(($nn + 1))
    jj=$(($ii + 9999))
    echo $nn $ii $jj
    if [[ $nn -ne 11 ]]
    then
        sed -e 's/FIRST/'$ii'/g' -e 's/LAST/'$jj'/g' -e 's/NN/'$nn'/g' template-trjcat > load.tcl
        vmd -dispdev text -e load.tcl
        rm -f load.tcl
    else
        echo no
    fi 
done
