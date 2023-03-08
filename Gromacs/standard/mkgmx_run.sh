#!/bin/bash

GPUID="$1"
GMX="gmx"
MDRUN="gmx mdrun -ntmpi 1 -ntomp 8 -gpu_id $GPUID"
MDRUN_GPU="gmx mdrun -ntmpi 1 -ntomp 8 -gpu_id $GPUID -pme gpu -nb gpu -bonded gpu -update gpu"
[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [GPU ID]; Suggested: 0"; exit 1; }

[ ! -d pdb2gmx ] && echo -e "mkgmx> Directory pdb2gmx not found!" && exit 1 

files=("pdb2gmx/ionized.gro" "pdb2gmx/topol.top" "pdb2gmx/index.ndx" )
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "mkgmx> $file not found!"
        exit 1
    fi
done

mkdir -p output
mini=true
double=false
heat=true
eq_npt=true
md1=true
md2=false

INITIAL_PDB=pdb2gmx/ionized.gro
NDX=pdb2gmx/index.ndx
TOP=pdb2gmx/topol.top

if $mini; then
    prefix=step1_mini
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f mdp/step1_mini.mdp -o $TPR -c $INITIAL_PDB -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -s $TPR -deffnm output/${prefix}
fi
    
if $double; then
    previous=step1_mini
    prefix=step1_mini_double
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f mdp/step1_mini_double.mdp -o $TPR -c output/${previous}.gro -r $INITIAL_PDB -n $NDX -p $TOP
    gmx_d mdrun -ntmpi 1 -ntomp 12 -s $TPR -v -deffnm output/${prefix}
fi

if $heat; then
    previous=step1_mini
    prefix=step3_annealing
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f mdp/step3_annealing.mdp -o $TPR -c output/${previous}.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -s $TPR -deffnm output/${prefix}
fi

if $eq_npt; then
    previous=step3_annealing
    prefix=step4_eq_npt
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f mdp/step4_eq_npt.mdp -o $TPR -c output/${previous}.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -s $TPR -deffnm output/${prefix}
fi

if $md1; then
    previous=step4_eq_npt
    prefix=md
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f mdp/step5_md.mdp -o $TPR -c output/${previous}.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN_GPU -v -s $TPR -deffnm output/${prefix}
fi

if $md2; then
    previous=md
    prefix=t1
    TPR=md.tpr
    $GMX grompp -f step5_md.mdp -o $TPR -t output/${previous}.cpt -c $INITIAL_PDB -r $INITIAL_PDB -n $NDX -p $TOP
    # gmx convert-tpr -s ${previous}.tpr -o $prefix.tpr -extend 10000
    # $MDRUN_GPU -v -s $TPR -cpi output/${prefix}.cpt -deffnm output/${prefix} -nsteps -1
    # $MDRUN_GPU -v -s $TPR -deffnm output/${prefix} -nsteps -1
fi
