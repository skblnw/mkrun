#!/bin/bash

MDRUN_GPU="gmx mdrun -ntmpi 1 -ntomp 16 -gpu_id 0 -pme gpu -nb gpu -bonded gpu -update gpu"

previous=t1
prefix=t1
TPR=md.tpr
# $GMX grompp -f step4_md.mdp -o output/step4_md.tpr -t output/$prefix1.cpt -n $NDX -p $TOP
# gmx convert-tpr -s $prefix1.tpr -o $prefix2.tpr -extend 10000
$MDRUN_GPU -v -s $TPR -cpi $previous.cpt -deffnm output/$prefix -nsteps -1
