#!/bin/bash

#gmx solvate -o solvate.pdb -box 10 10 10 -maxsol 32000
gmx solvate -o solvate.pdb -box 17 17 17 -maxsol 160000

gmx pdb2gmx -f solvate.pdb -o tmp.pdb <<EOF
8
1
EOF

gmx grompp -f md.mdp -o md.tpr -c solvate.pdb -p topol.top

rm -f solvate.pdb topol.top tmp.pdb
