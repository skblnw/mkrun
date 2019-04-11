#!/bin/bash

mkdir -p output
mini1=true
mini2=true
heat=true
eq_npt=true
md1=false
md2=false

if $mini1; then
    prefix=step1_mini_crystal
    gmx grompp -f $prefix.mdp -o $prefix.tpr -c pdb2gmx/ionized.pdb -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top -maxwarn 1
    gmx mdrun -pin on -ntmpi 1 -ntomp 8 -v -s $prefix.tpr -deffnm output/$prefix
fi
    
if $mini2; then
    prefix1=step1_mini_crystal
    prefix2=step2_mini
    gmx_d grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top  -maxwarn 1
    gmx_d mdrun -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $heat; then
    prefix1=step2_mini
    prefix2=step3_annealing
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top  -maxwarn 1
    gmx mdrun -pin on  -ntmpi 1 -ntomp 8 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $eq_npt; then
    prefix1=step3_annealing
    prefix2=step4_eq_npt
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top  -maxwarn 1
    gmx mdrun -pin on -ntmpi 1 -ntomp 8 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $md1; then
    prefix1=step4_eq_npt
    prefix2=step5_md-1
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top  -maxwarn 1
    gmx mdrun -pin on -ntmpi 1 -ntomp 8 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $md2; then
    prefix1=step5_md-1
    prefix2=step5_md-2
#    gmx convert-tpr -s $prefix1.tpr -o $prefix2.tpr -extend 10000
    gmx mdrun -pin on -ntmpi 1 -ntomp 8 -v -noappend -s step5_md-1.tpr -cpi output/$prefix1.cpt -deffnm output/$prefix2 -nsteps -1
fi
