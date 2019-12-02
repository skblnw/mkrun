#!/bin/bash

# Short equilibration
gmx grompp -f npt_umbrella.mdp -c gro/confXXX.gro -p pdb2gmx/topol.top -n pdb2gmx/index.ndx -o output/nptXXX.tpr -maxwarn 1
gmx mdrun -ntomp 8 -pin on -v -deffnm output/nptXXX

# Umbrella run
gmx grompp -f md_umbrella.mdp -c output/nptXXX.gro -t output/nptXXX.cpt -p pdb2gmx/topol.top -n pdb2gmx/index.ndx -o output/umbrellaXXX.tpr
gmx mdrun -ntomp 8 -pin on -v -deffnm output/umbrellaXXX 

