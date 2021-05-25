#!/bin/bash

gmx solvate -o solvate.pdb -box 10 10 10 -maxsol 32000
#gmx solvate -o solvate.pdb -box 12 12 12 -maxsol 55000
#gmx solvate -o solvate.pdb -box 17 17 17 -maxsol 160000

gmx pdb2gmx -f solvate.pdb <<EOF
8
1
EOF

gmx grompp -f mini.mdp -o run_mini.tpr -c solvate.pdb -p topol.top
gmx mdrun -v -deffnm run_mini
gmx grompp -f md.mdp -o run_md.tpr -c run_mini.gro -p topol.top

rm -f solvate.pdb mdout.mdp run_mini.*
