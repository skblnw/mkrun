#!/bin/bash

rm -f ionized*

cat > tcl <<'EOF'
package require Orient
set filename ../prot
mol new $filename.pdb
set sel [atomselect top "all"]
set I [Orient::calc_principalaxes $sel]
set A [Orient::orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [Orient::calc_principalaxes $sel]
set A [Orient::orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
$sel writepdb tmp.pdb

package require solvate
#solvate $filename.psf tmp.pdb -b 2.4 -minmax {{-77 -76 -35} {77 76 130}} -o solvated
solvate $filename.psf tmp.pdb -b 2.4 -t 10 -o solvated
package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -cation POT -o ionized
quit
EOF

vmd -dispdev text -e tcl
rm -f tcl tmp.pdb solvated.p*
