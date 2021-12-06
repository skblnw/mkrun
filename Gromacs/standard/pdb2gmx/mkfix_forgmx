#!/bin/bash
# 
# ZN2 for ZN
# set sel [atomselect top "resname ZN"]
# \$sel set resname ZN2
# \$sel set chain Z
# \$sel writepdb chains/ZN.pdb
#
# CLA for CL
# set sel [atomselect top "resname CL"]
# \$sel set resname CLA
# \$sel set name CLA
# \$sel set chain Y
# \$sel writepdb chains/CL.pdb
#
VMD="/opt/vmd/1.9.3/vmd"

PDB="../../../../raw/6m0j.pdb"

rm -r chains mol.p*
mkdir chains
cat > tcl <<EOF
mol new $PDB
foreach ii {A E} {
  set sel [atomselect top "protein and chain \$ii"]
  \$sel writepdb chains/\$ii.pdb
}
set sel [atomselect top "resname ZN"]
\$sel set resname ZN2
\$sel set chain Z
\$sel writepdb chains/ZN.pdb
set sel [atomselect top "resname CL"]
\$sel set resname CLA
\$sel set name CLA
\$sel set chain Y
\$sel writepdb chains/CL.pdb
set sel [atomselect top "resname HOH"]
\$sel set chain X
\$sel writepdb chains/HOH.pdb
quit
EOF
$VMD -dispdev text -e tcl 

for ii in A E ZN CL HOH; do
  grep "^ATOM" chains/$ii.pdb >> mol.pdb
done

cat > tcl <<'EOF'
mol new mol.pdb
set sel [atomselect top "resname HSD HSE HSP"]; $sel set resname HIS
set sel [atomselect top "resname ILE and name CD"]; $sel set name CD1
set sel [atomselect top "not hydrogen"]
$sel writepdb mol.pdb
quit
EOF
$VMD -dispdev text -e tcl

rm tcl
