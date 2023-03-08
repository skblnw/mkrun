#!/bin/bash

MOL="$1"
NOMOVE=true
[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [MOLECULE]; Suggested: conf.gro"; exit 1; }

if [ ! -f $MOL ]; then
    echo "mkgmx> $MOL \nStarting system not found!"
    exit 1
fi

rm editconf.gro solvate.gro ions.tpr ionized.gro

cat > ions.mdp << EOF
integrator  = md
dt          = 0.001
nsteps      = 50000
nstlist         = 1
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF

#gmx pdb2gmx -f mol.pdb

# Create a cubic box of 2 nm on each side and place the system to the center and align its principle axes to the reference axes
# -princ usually helps to reduce your system size
echo 0 | gmx editconf -f $MOL -o editconf.gro -princ -d 1.0 -bt cubic
#echo 0 | gmx editconf -f $MOL -o editconf.gro -princ -center 3.7 3.4 2.4 -box 8.0 12.2 5.2

# Solvate the box
gmx solvate -cp editconf.gro -o solvate.gro -p topol.top
# Add ions to make it neutral and of 0.15 M NaCl
# If you want KCl, add -pname K
gmx grompp -f ions.mdp -c solvate.gro -o ions.tpr -p topol.top -maxwarn 1 > LOG_grompp 2>&1
gmx genion -s ions.tpr -o ionized.gro -conc 0.15 -neutral -p topol.top
echo q | gmx make_ndx -f ionized.gro 

if $NOMOVE; then
# Add NOMOVE
cp topol_Protein_chain_A.itp topol_Protein_chain_A.itp.BAK
cat >> topol_Protein_chain_A.itp <<EOF

#ifdef POSRES_NOMOVE
#include "posre_Protein_chain_A_NOMOVE.itp"
#endif
EOF

cat > posre_Protein_chain_A_NOMOVE.itp <<EOF
[ position_restraints ]
     3199     1  1000  1000  1000
EOF
fi

rm \#* editconf.gro solvate.gro ions.tpr ions.mdp mdout.mdp
