#!/bin/bash

jj=0
for ii in $(seq 50.0 1.5 62.0)
do
    echo $jj $ii 
    mkdir -p $ii 
    cp ../pulling/50.0.pdb $ii/mol.pdb
    cp fix.tcl *.in COM_dist.RST $ii
    cd $ii; vmd -dispdev text -e fix.tcl; sed -i "s/DISTHERE/${ii}/g" COM_dist.RST; tleap -f ../build.leap; cd ..

    jj=$((jj+1))
done