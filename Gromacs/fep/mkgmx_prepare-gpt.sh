#!/bin/bash

#########################################
# Description: Give you a peptide_mutate.pdb in which 1. mutated residues were renamed 2. missing atoms were added (if necessary).
# This runs before gentop.
# Authors:
# Kevin C. Chan (work@skblnw.com) 2022.10.25
#########################################

# Set variables
CHAIN_DIR="./split/chains"
LOG_GROMPP="LOG_grompp"
IONS_MDP="ions.mdp"

# Create symlinks
ln -fs $CHAIN_DIR chains

# Define function to cleanup temporary files
cleanup() {
  rm -f \#* editconf.pdb solvate.pdb ions.mdp ions.tpr index.ndx mdout.mdp LOG_grompp
  rm -f conf.pdb free.pdb ionized.pdb mol.pdb peptide_mutate.pdb X.pdb
}

# Remove previous topol.top file
rm -f topol.top

# Run gmx to get conf.gro
gmx pdb2gmx -f chains/chain_C.pdb -ff charmm36-feb2021 -water tip3p

# Mutate conf.gro using pmx
pmx mutate -f conf.gro -ff charmm36m-mut -o peptide_mutate.pdb --script mut

# Divide (and fix the names of) complex PDB into multiple fragments using segname
cat > tcl <<'EOF'
mol new peptide_mutate.pdb
set sel [atomselect top "all"]
$sel set chain X
$sel writepdb X.pdb
quit
EOF
vmd -dispdev text -e tcl
rm tcl

# Get free.pdb
rm -f mol.pdb
grep "^ATOM" X.pdb > mol.pdb
gmx pdb2gmx -f mol.pdb -ff charmm36m-mut -water tip3p -o free.pdb

# Get conf.pdb (complex) and topol.top (complex)
for ii in A B; do
    grep "^ATOM" chains/chain_$ii.pdb >> mol.pdb
done
gmx pdb2gmx -f mol.pdb -ff charmm36m-mut -water tip3p -o conf.pdb

# Generate hybrid topology using pmx
pmx gentop -p topol.top -ff charmm36m-mut
sed 's/^Protein_chain_[AB].*//g' pmxtop.top > pmxtop_free.top

# Create ions.mdp file
cat > $IONS_MDP << EOF
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

# Remove and create directories
rm -rf free complex
mkdir -p free/pdb2gmx complex/pdb2gmx

# Function for creating box, solvating, and adding ions
prepare_system() {
  local system_name=$1
  local input_pdb=$2
  local output_top=$3

  # Create the box
  rm -f editconf.pdb solvate.pdb ions.tpr index.ndx
  echo 0 | gmx editconf -f $input_pdb -o editconf.pdb -princ -d 1.0 -bt cubic

  # Solvate the box
  gmx solvate -cp editconf.pdb -o solvate.pdb -p $output_top

  # Add ions to make it neutral and of 0.15 M NaCl
  # If you want KCl, add -pname K
  gmx grompp -f $IONS_MDP -c solvate.pdb -o ions.tpr -p $output_top -maxwarn 1 > $LOG_GROMPP 2>&1
  echo 13 | gmx genion -s ions.tpr -o ionized.pdb -conc 0.15 -neutral -p $output_top

  # Create index file
  echo q | gmx make_ndx -f ionized.pdb

  # Copy files to the respective directory
  cp ionized.pdb index.ndx pmx_* posre_* $system_name/pdb2gmx
  cp $output_top $system_name/pdb2gmx/topol.top
}

# Prepare complex system
prepare_system "complex" "conf.pdb" "pmxtop.top"

# Prepare free system
# Get the box size first
gmx editconf -f ionized.pdb -o tmp.gro
bs=$(tail -1 tmp.gro)
rm ionized.pdb tmp.gro

prepare_system "free" "free.pdb" "pmxtop_free.top" "$bs"

# Cleanup temporary files
cleanup
