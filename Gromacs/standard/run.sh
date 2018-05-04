#!/bin/bash

if true; then
    prefix=step1_mini_crystal
    gmx grompp -f $prefix.mdp -o $prefix.tpr -c pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -ntomp 8 -v -s $prefix.tpr -deffnm output/$prefix
fi
    
if true; then
    prefix1=step1_mini_crystal
    prefix2=step2_mini
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -ntomp 8 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if true; then
    prefix1=step2_mini
    prefix2=step3_annealing
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -ntomp 8 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if true; then
    prefix1=step3_annealing
    prefix2=step4_eq_npt
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -ntomp 8 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if true; then
    prefix1=step4_eq_npt
    prefix2=step5_md
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx_cuda mdrun -ntomp 8 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if false; then
    prefix1=step5_md
    prefix2=step5_md-2
    gmx convert-tpr -s $prefix1.tpr -o $prefix2.tpr -extend 1000
    gmx mdrun -ntomp 8 -v -s $prefix2.tpr -cpi $prefix1.cpt -deffnm output/$prefix2
fi