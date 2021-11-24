#!/bin/bash

rm -f prot*
DIR="../../7rtr_charmmgui/"

cat > tcl <<EOF
package require topotools
mol new $DIR/step1_pdbreader.psf
mol addfile $DIR/step1_pdbreader.pdb
set sel [atomselect top "segname PROA PROB PROC"]
set mol [::TopoTools::selections2mol "\$sel"]
animate write psf prot.psf \$mol
animate write pdb prot.pdb \$mol
quit
EOF

vmd -dispdev text -e tcl 
rm tcl