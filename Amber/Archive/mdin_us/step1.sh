#!/bin/bash

jj=0
for ii in $(seq 32.0 1.5 60.5)
do
    ipdb=$(($jj * 7))
    echo $jj $ii $ipdb
    rm -rf $ii
    mkdir -p $ii 
    cp ../pulling/${ipdb}.pdb $ii/mol.pdb
    cd $ii; vmd -dispdev text mol.pdb -e ../fix.tcl; tleap -f ../build.leap; cd ..

    jj=$(($jj+1))
done
