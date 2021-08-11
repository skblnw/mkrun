#!/bin/bash

NAMD="namd3 +p1 +devices 1"

for ii in $(seq 1 5); do
    rsync -avh eq/fep.namd eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc t$ii
    $NAMD t$ii/fep.namd >& t$ii/LOG_fep
done
