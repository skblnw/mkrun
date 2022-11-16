#!/bin/bash
GMX="gmx"
MDRUN="gmx mdrun -ntmpi 1 -ntomp 8 -gpu_id 0"
MDRUN_GPU="gmx mdrun -ntmpi 1 -ntomp 8 -gpu_id 0 -pme gpu -nb gpu -bonded gpu -update gpu"

mkdir -p output
mini=true
double=false
heat=true
eq_npt=true
md1=true
md2=false

[ ! -d pdb2gmx ] && echo -e "mkgmx> Directory pdb2gmx not found!" && exit 1 
[ ! -f $NDX ] && echo -e "mkgmx> Index file not found!" && exit 1 

INITIAL_PDB=pdb2gmx/ionized.pdb
NDX=pdb2gmx/index.ndx
TOP=pdb2gmx/topol.top

if $mini; then
    prefix=step1_mini
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f step1_mini.mdp -o $TPR -c $INITIAL_PDB -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -s $TPR -deffnm output/${prefix}
fi
    
if $double; then
    previous=step1_mini
    prefix=step1_mini_double
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f step1_mini_double.mdp -o $TPR -c output/${previous}.gro -r $INITIAL_PDB -n $NDX -p $TOP
    gmx_d mdrun -ntmpi 1 -ntomp 12 -s $TPR -v -deffnm output/${prefix}
fi

if $heat; then
    previous=step1_mini
    prefix=step3_annealing
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f step3_annealing.mdp -o $TPR -c output/${previous}.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -s $TPR -deffnm output/${prefix}
fi

if $eq_npt; then
    previous=step3_annealing
    prefix=step4_eq_npt
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f step4_eq_npt.mdp -o $TPR -c output/${previous}.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -s $TPR -deffnm output/${prefix}
fi

if $md1; then
    previous=step4_eq_npt
    prefix=md
    TPR=${prefix}.tpr
    rm -f $TPR
    $GMX grompp -f step5_md.mdp -o $TPR -c output/${previous}.gro -n $NDX -p $TOP
    $MDRUN_GPU -v -s $TPR -deffnm output/${prefix}
fi

if $md2; then
    previous=md
    prefix=t1
    TPR=md.tpr
    # $GMX grompp -f step4_md.mdp -o output/step4_md.tpr -t output/${previous}.cpt -n $NDX -p $TOP
    # gmx convert-tpr -s ${previous}.tpr -o $prefix.tpr -extend 10000
    $MDRUN_GPU -v -s $TPR -cpi output/${prefix}.cpt -deffnm output/${prefix} -nsteps -1
fi
