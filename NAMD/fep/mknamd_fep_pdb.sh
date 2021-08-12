#!/bin/bash

cat > tcl <<'EOF'
mol new pdb2namd/vmd_solvate/ionized.pdb
for {set ii 1} {$ii < 11} {incr ii} {
    mol addfile t$ii/alchemy.coor
    set sel [atomselect top all]
    $sel writepdb $ii.pdb
}
quit
EOF

vmd -dispdev text -e tcl
rm tcl