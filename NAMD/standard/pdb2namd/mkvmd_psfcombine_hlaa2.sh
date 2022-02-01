#!/bin/bash
VMD="/opt/vmd/1.9.3/vmd"

rm -rf chains peptide.p* prot* 
mkdir chains
DIR1="../../raw/7kgq.charmm-gui-4304186123"
DIR2="../../../1oga/raw/1oga.tcr.charmm-gui-4305997050"

cat > tcl <<EOF
package require topotools
mol new $DIR1/step1_pdbreader.psf
mol addfile $DIR1/step1_pdbreader.aligned.pdb
mol new $DIR2/step1_pdbreader.psf
mol addfile $DIR2/step1_pdbreader.aligned.pdb
set sel1 [atomselect 0 "segname PROA PROB PROC"]
set sel2 [atomselect 1 "segname PROD and resid < 112 or segname PROE and resid < 116"]
set mol [::TopoTools::selections2mol "\$sel1 \$sel2"]
animate write psf prot.psf \$mol
animate write pdb prot.pdb \$mol
quit
EOF
$VMD -dispdev text -e tcl 
rm tcl
