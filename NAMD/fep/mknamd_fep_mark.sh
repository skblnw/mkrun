#!/bin/bash

cat > tcl <<EOF
mol new pdb2namd/vmd_solvate/ionized.pdb
set sel [atomselect top all]
\$sel set beta 0
set sel [atomselect top "protein and not {backbone or name HA HB CB} and name \".*A\""]
\$sel set beta -1
set sel [atomselect top "protein and not {backbone or name HA HB CB} and name \".*B\""]
\$sel set beta 1
set sel [atomselect top all]
\$sel writepdb ionized.fep
quit
EOF

vmd -dispdev text -e tcl
rm tcl