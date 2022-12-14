#!/bin/bash

## Get the peptide using segname
mkdir chains
cat > tcl <<'EOF'
mol new label.pdb
set sel [atomselect top "resname HSD HSE HSP"]; $sel set resname HIS
set sel [atomselect top "resname ILE and name CD"]; $sel set name CD1
set sel [atomselect top "name OT1"]; $sel set name O
set sel [atomselect top "name OT2"]; $sel set name OXT
set sel [atomselect top "noh chain A"]; $sel writepdb chains/A.pdb
set sel [atomselect top "noh chain B"]; $sel writepdb chains/B.pdb
set sel [atomselect top "noh chain X"]; $sel writepdb chains/C.pdb
quit
EOF
vmd -dispdev text -e tcl
rm tcl
