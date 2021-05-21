#!/bin/bash
GMX="gmx"
MDRUN="gmx mdrun -ntmpi 1 -ntomp 16 -gpu_id 1"

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
    TPR=output/$prefix.tpr
    rm -f $TPR
    $GMX grompp -f $prefix.mdp -o $TPR -c $INITIAL_PDB -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -deffnm output/$prefix
fi
    
if $double; then
    prefix1=step1_mini
    prefix2=step1_mini_double
    TPR=output/$prefix2.tpr
    rm -f $TPR
    $GMX grompp -f $prefix2.mdp -o $TPR -c output/$prefix1.gro -r $INITIAL_PDB -n $NDX -p $TOP
    gmx_d mdrun -ntmpi 1 -ntomp 12 -v -deffnm output/$prefix2
fi

if $heat; then
    prefix1=step1_mini
    prefix2=step3_annealing
    TPR=output/$prefix2.tpr
    rm -f $TPR
    $GMX grompp -f $prefix2.mdp -o $TPR -c output/$prefix1.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -deffnm output/$prefix2
fi

if $eq_npt; then
    prefix1=step3_annealing
    prefix2=step4_eq_npt
    TPR=output/$prefix2.tpr
    rm -f $TPR
    $GMX grompp -f $prefix2.mdp -o $TPR -c output/$prefix1.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -deffnm output/$prefix2
fi

if $md1; then
    prefix1=step4_eq_npt
    prefix2=step5_md
    TPR=output/$prefix2.tpr
    rm -f $TPR
    $GMX grompp -f $prefix2.mdp -o $TPR -c output/$prefix1.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -deffnm output/$prefix2
fi

if $md2; then
    prefix=output/t1
#    gmx convert-tpr -s $prefix1.tpr -o $prefix2.tpr -extend 10000
    $MDRUN -v -noappend -s output/step5_md.tpr -cpi output/step5_md.cpt -deffnm $prefix -nsteps -1
fi
