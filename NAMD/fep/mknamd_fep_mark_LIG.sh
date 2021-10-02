#!/bin/bash

cat > tcl <<EOF
mol new pdb2namd/vmd_solvate/ionized.pdb
set sel [atomselect top all]
\$sel set beta 0
set sel [atomselect top "resname LIG and name C4A C6A C8A C9A H7A H8A"]
\$sel set beta -1
set sel [atomselect top "resname LIG and name C4B C6B C8B C9B H7B H8B"]
\$sel set beta 1
set sel [atomselect top all]
\$sel writepdb ionized.fep
quit
EOF

vmd -dispdev text -e tcl
rm tcl