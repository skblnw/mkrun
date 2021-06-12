#!/bin/bash

NAMD="namd3 +p16 +devices 0"

for ii in $(seq 1 5); do
    rsync -avh eq/ t$ii --exclude="*.BAK" --exclude="*.old" --exclude="*.fepout" --exclude="LOG_eq"
    cd t$ii
    $NAMD fep.soft.namd >& LOG_fep
    cd ..
done