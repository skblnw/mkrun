#!/bin/bash

rm ionized*

cat > tcl <<'EOF'
set filename ../prot

package require Orient
mol new $filename.pdb
set sel [atomselect top "all"]
set I [Orient::calc_principalaxes $sel]
set A [Orient::orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [Orient::calc_principalaxes $sel]
set A [Orient::orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
$sel moveby [vecinvert [measure center $sel]]
$sel writepdb tmp.pdb

package require solvate
solvate $filename.psf tmp.pdb -t 12 -o solvated
# solvate $filename.psf tmp.pdb -minmax {{-48 -48 -70} {48 48 70}} -o solvated
# solvate $filename.psf tmp.pdb -minmax {{-48 -48 -48} {48 48 48}} -o solvated
package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -o ionized

mol new ionized.pdb type pdb waitfor all
set all [atomselect top "all"]
set fout [open "cell_size.str" w]
set all [atomselect top water] 
set minmax [measure minmax $all] 
set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]] 
puts $fout "cellBasisVector1 [lindex $vec 0] 0 0" 
puts $fout "cellBasisVector2 0 [lindex $vec 1] 0" 
puts $fout "cellBasisVector3 0 0 [lindex $vec 2]" 
set center [measure center $all] 
puts $fout "cellOrigin $center" 
close $fout
quit
EOF

vmd -dispdev text -e tcl
rm tcl tmp.pdb solvated.p*
