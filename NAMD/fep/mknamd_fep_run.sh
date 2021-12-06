#!/bin/bash

NAMD="namd3 +p1 +devices 1"

for ii in $(seq 1 5); do
    rsync -avh eq/fep.namd eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc trial$ii
    $NAMD trial$ii/fep.namd >& trial$ii/LOG_fep
done
