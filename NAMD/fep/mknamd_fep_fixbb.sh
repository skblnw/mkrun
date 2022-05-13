#!/bin/bash

cat > tcl <<EOF
mol new pdb2namd/vmd_solvate/ionized.pdb type pdb waitfor all
set all [atomselect top "all"]
\$all set beta 0
set sel [atomselect top "segname MUT and name C N O"]
\$sel get resid
\$sel get resname
\$sel get name
\$sel set beta 1
\$all writepdb cons.fep
quit
EOF

vmd -dispdev text -e tcl >& /dev/null; rm tcl
