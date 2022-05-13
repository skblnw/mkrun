#!/bin/bash

NAMD="namd3"
DEV=1

cd eq
$NAMD +p8 +devices $DEV fep.eq.namd >& LOG_eq
cd ..

for ii in $(seq 1 1); do
    rsync -avh eq/fep.namd eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc trial$ii
    $NAMD +p1 +devices $DEV trial$ii/fep.namd >& trial$ii/LOG_fep
done
