#!/bin/bash

echo 0 | gmx trjconv -f ed.trr -s ../ed.tpr -split 20

for ii in $(seq 0 10); do
    jj=$(($ii + 1))
    echo 0 | gmx trjconv -f trajout$ii.xtc -s ../ed.tpr -skip 20 -o $jj.gro
    rm -f trajout$ii.xtc
done