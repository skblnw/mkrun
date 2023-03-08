#!/bin/bash

GPUID="$1"
PREFIX="$2"
PREVIOIUS="$3"
MDRUN_GPU="gmx mdrun -ntmpi 1 -ntomp 16 -gpu_id $GPUID -pme gpu -nb gpu -bonded gpu -update gpu"
[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [GPU ID] [PEFIX] [PREVIOIUS]; Suggested: 0 t1; [PREVIOIUS] can be empty"; exit 1; }

[ ! -f "md.tpr" ] && echo -e "mkgmx> md.tpr not found!" && exit 1 

prefix=$PREFIX
tpr=md.tpr

if [ -z "$PREVIOIUS" ]; then
    [ -f "${PREFIX}.cpt" ] && echo -e "mkgmx> ${PREFIX}.cpt found but you wanna start a new one." && exit 1 
    echo "mkvmd> Start a new simulation"
    $MDRUN_GPU -v -s ${tpr} -deffnm ${prefix} -nsteps -1
else
    [ ! -f "${PREVIOIUS}.cpt" ] && echo -e "mkgmx> ${PREVIOIUS}.cpt not found!" && exit 1 
    cp ${PREVIOIUS}.cpt ${PREVIOIUS}.cpt.BAK
    $MDRUN_GPU -v -s ${tpr} -cpi ${PREVIOIUS}.cpt -deffnm ${prefix} -nsteps -1
fi

# gmx grompp -o md.tpr -f step7_production_vrescale.mdp -c step7_1.gro -t step7_1.cpt -n index.ndx -p topol.top
# gmx convert-tpr -s $prefix1.tpr -o $prefix2.tpr -extend 10000
