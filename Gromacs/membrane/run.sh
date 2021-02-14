#!/bin/bash

mkdir -p output
em1=false
em2=false
eq_nvt1=false
eq_npt1=false
eq_npt2=true
eq_npt3=true
eq_npt4=true
eq_npt5=true
eq_npt6=true
md1=false
md2=false

if $em1; then
    prefix=step1.0_em
    rm -f $prefix.tpr
    gmx grompp -f $prefix.mdp -o $prefix.tpr -c pdb2gmx/ionized.pdb -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -pin on -ntmpi 1 -ntomp 12 -v -s $prefix.tpr -deffnm output/$prefix
fi
    
if $em2; then
    prefix1=step1.0_em
    prefix2=step1.1_em
    rm -f $prefix2.tpr
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -ntmpi 1 -ntomp 12 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $eq_nvt1; then
    prefix1=step1.1_em
    prefix2=step2.0_eq_nvt
    rm -f $prefix2.tpr
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun  -ntmpi 1 -ntomp 12 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $eq_npt1; then
    prefix1=step2.0_eq_nvt
    prefix2=step3.0_eq_npt
    rm -f $prefix2.tpr
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -pin on -ntmpi 1 -ntomp 12 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $eq_npt2; then
    prefix1=step3.0_eq_npt
    prefix2=step3.1_eq_npt
    rm -f $prefix2.tpr
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -pin on -ntmpi 1 -ntomp 12 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $eq_npt3; then
    prefix1=step3.1_eq_npt
    prefix2=step3.2_eq_npt
    rm -f $prefix2.tpr
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -pin on -ntmpi 1 -ntomp 12 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $eq_npt4; then
    prefix1=step3.2_eq_npt
    prefix2=step3.3_eq_npt
    rm -f $prefix2.tpr
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -pin on -ntmpi 1 -ntomp 12 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $eq_npt5; then
    prefix1=step3.3_eq_npt
    prefix2=step3.4_eq_npt
    rm -f $prefix2.tpr
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -pin on -ntmpi 1 -ntomp 12 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $eq_npt6; then
    prefix1=step3.4_eq_npt
    prefix2=step3.5_eq_npt
    rm -f $prefix2.tpr
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -pin on -ntmpi 1 -ntomp 12 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $md1; then
    prefix1=step3.5_eq_npt
    prefix2=step4_md-1
    rm -f $prefix2.tpr
    gmx grompp -f $prefix2.mdp -o $prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx mdrun -pin on -ntmpi 1 -ntomp 12 -v -s $prefix2.tpr -deffnm output/$prefix2
fi

if $md2; then
    prefix1=step5_md-1
    prefix2=step5_md-2
#    gmx convert-tpr -s $prefix1.tpr -o $prefix2.tpr -extend 10000
    gmx mdrun -pin on -ntmpi 1 -ntomp 12 -v -noappend -s step5_md-1.tpr -cpi output/$prefix1.cpt -deffnm output/$prefix2 -nsteps -1
fi
