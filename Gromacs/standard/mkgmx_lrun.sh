#!/bin/bash

# Check if the number of arguments is zero
if [ $# -eq 0 ]; then 
    echo "mkvmd> Usage: $0 [GPU ID] [PREFIX] [PREVIOUS]"
    echo "mkvmd> Suggested: 0 t1; [PREVIOUS] can be empty"
    echo "mkvmd> OR Suggested: 0 t1 t1; [PREVIOUS] only contains the prefix"
    exit 1
fi

# Define parameters
GPUID="$1"
PREFIX="$2"
PREVIOUS="$3"
TPR_FILE="md.tpr"
MDRUN_GPU="gmx mdrun -ntmpi 1 -ntomp 16 -gpu_id $GPUID -pme gpu -nb gpu -bonded gpu -update gpu"

# Check if the TPR file exists
[ ! -f "$TPR_FILE" ] && echo "mkvmd> $TPR_FILE not found!" && exit 1 

# Start a new simulation or continue from previous checkpoint
if [ -z "$PREVIOUS" ]; then
    [ -f "${PREFIX}.cpt" ] && echo "mkvmd> ${PREFIX}.cpt found but you want to start a new one." && exit 1 
    echo "mkvmd> Starting a new simulation"
    $MDRUN_GPU -v -s $TPR_FILE -deffnm $PREFIX -nsteps -1
else
    [ ! -f "${PREVIOUS}.cpt" ] && echo "mkvmd> ${PREVIOUS}.cpt not found!" && exit 1 
    cp ${PREVIOUS}.cpt ${PREVIOUS}.cpt.BAK
    echo "mkvmd> Continuing simulation from ${PREVIOUS}.cpt"
    $MDRUN_GPU -v -s $TPR_FILE -cpi ${PREVIOUS}.cpt -deffnm $PREFIX -nsteps -1
fi

# gmx grompp -o md.tpr -f step7_production_vrescale.mdp -c step7_1.gro -t step7_1.cpt -n index.ndx -p topol.top
# gmx convert-tpr -s $prefix1.tpr -o $prefix2.tpr -extend 10000
