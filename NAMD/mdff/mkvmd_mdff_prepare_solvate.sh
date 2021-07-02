#!/bin/bash

# Must have:
# 1. prot.psf
# 2. prot.pdb
# 3. map.mrc

source ~/.zshrc

PREFIX=ionized
PSF=${PREFIX}.psf
PDB=${PREFIX}.pdb
MAP=${PREFIX}.mrc

rm -f ${PREFIX}_grid.pdb ${PREFIX}_ssrestraints.txt ${PREFIX}_cispeptide.txt ${PREFIX}_chirality.txt ${PREFIX}_grid.dx 

cat > tmp.tcl << EOF

package require mdff
package require ssrestraints
package require cispeptide
package require chirality

mdff griddx -i $MAP -o tmp.dx
mdff gridpdb -psf $PSF -pdb $PDB -o ${PREFIX}_grid.pdb
ssrestraints -psf $PSF -pdb $PDB -o ${PREFIX}_ssrestraints.txt -hbonds

mol new $PSF
mol addfile $PDB
set sel [atomselect top all]
cispeptide restrain -o ${PREFIX}_cispeptide.txt
chirality restrain -o ${PREFIX}_chirality.txt

set min [vecsub [lindex [measure minmax \$sel] 0] {1 1 1}]
set max [vecadd [lindex [measure minmax \$sel] 1] {1 1 1}]
set amt [list [expr int([lindex \$min 0])] [expr int([lindex \$min 1])] [expr int([lindex \$min 2])] [expr int([lindex \$max 0])] [expr int([lindex \$max 1])] [expr int([lindex \$max 2])]]
voltool crop -amt \$amt -i tmp.dx -o ${PREFIX}_grid.dx

mdff setup -o mdff -psf $PSF -pdb $PDB -griddx ${PREFIX}_grid.dx -gridpdb ${PREFIX}_grid.pdb -extrab {${PREFIX}_ssrestraints.txt ${PREFIX}_cispeptide.txt ${PREFIX}_chirality.txt} \\
-gscale 0.1 -minsteps 5000 -numsteps 10000 -pbc

quit
EOF

vmd -dispdev text -e tmp.tcl 

cp mdff-step1.namd mdff-step1.namd.BAK

sed -e "s/parameters par_all36_prot.prm//g" \
-e "s/paraTypeCharmm on/paraTypeCharmm on\nparameters par_all36m_prot.prm\nmergeCrossterms yes\nparameters par_all36_lipid.prm\nparameters par_all36_carb.prm\nparameters par_all36_na.prm\nparameters par_all36_cgenff.prm\nparameters toppar_water_ions_namd.str\n\nmgridforcechecksize 0 off/g" \
mdff-step1.namd.BAK > mdff-step1.namd

rm -f tmp.tcl tmp.dx
