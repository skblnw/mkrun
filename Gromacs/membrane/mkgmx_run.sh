#!/bin/bash
GMX="gmx"
MDRUN="gmx mdrun -ntmpi 1 -ntomp 8 -gpu_id 1"
MDRUN_GPU="gmx mdrun -ntmpi 1 -ntomp 8 -gpu_id 1 -pme gpu -nb gpu -bonded gpu -update gpu"

mkdir -p output
mini=false
double=false
heat=false
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
    prefix2=step2.0_eq_nvt
    TPR=output/$prefix2.tpr
    rm -f $TPR
    $GMX grompp -f $prefix2.mdp -o $TPR -c output/$prefix1.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -deffnm output/$prefix2
fi

if $eq_npt; then
    for ii in $(seq 0 5); do
        if [[ $ii -eq 0 ]]; then
            prefix1=step2.0_eq_nvt
            prefix2=step3.${ii}_eq_npt
        else
            jj=$((ii-1))
            prefix1=step3.${jj}_eq_npt
            prefix2=step3.${ii}_eq_npt
        fi
        TPR=output/$prefix2.tpr
        rm -f $TPR
        $GMX grompp -f $prefix2.mdp -o $TPR -c output/$prefix1.gro -r $INITIAL_PDB -n $NDX -p $TOP
        $MDRUN -v -deffnm output/$prefix2
    done
fi

if $md1; then
    prefix1=step3.5_eq_npt
    prefix2=step4_md
    TPR=output/$prefix2.tpr
    rm -f $TPR
    $GMX grompp -f $prefix2.mdp -o $TPR -c output/$prefix1.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -deffnm output/$prefix2
fi

if $md2; then
    prefix=output/t1
#    gmx convert-tpr -s $prefix1.tpr -o $prefix2.tpr -extend 10000
    $MDRUN_GPU -v -noappend -s output/step4_md-cont.tpr -cpi output/step4_md.cpt -deffnm $prefix -nsteps -1
fi
