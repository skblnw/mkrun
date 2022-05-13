#!/bin/bash

GPUID=0
MDRUN_GPU="gmx mdrun -ntmpi 1 -ntomp 16 -gpu_id $GPUID -pme gpu -nb gpu -bonded gpu -update gpu"

previous=t1
prefix=t1
TPR=md.tpr
# gmx grompp -o md.tpr -f step7_production_vrescale.mdp -c step7_1.gro -t step7_1.cpt -n index.ndx -p topol.top
# gmx convert-tpr -s $prefix1.tpr -o $prefix2.tpr -extend 10000
$MDRUN_GPU -v -s $TPR -cpi $previous.cpt -deffnm $prefix -nsteps -1
