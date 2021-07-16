#!/bin/bash

MUT="452 478"

cat > tcl <<EOF
mol new pdb2namd/vmd_solvate/ionized.pdb
set sel [atomselect top all]
\$sel set beta 0
foreach ii {$MUT} {
    set sel [atomselect top "resid \$ii and not {backbone or name HA} and name \".*A\""]
    \$sel set beta -1
    set sel [atomselect top "resid \$ii and not {backbone or name HA} and name \".*B\""]
    \$sel set beta 1
}
set sel [atomselect top all]
\$sel writepdb ionized.fep
quit
EOF

vmd -dispdev text -e tcl
rm tcl