#!/bin/bash

NAMD="/opt/namd/NAMD_3.0alpha8_Linux-x86_64-multicore-CUDA/namd3 +p4 +devices 0"

for ii in $(seq 1 5); do
    rsync -avh eq/ t$ii --exclude="*.BAK" --exclude="*.old" --exclude="*.fepout" --exclude="LOG_eq"
    cd t$ii
    $NAMD fep.soft.namd >& LOG_fep
    cd ..
done