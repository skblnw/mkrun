#!/bin/bash

SEL_MHC="segname PROA"
SEL_PEP="segname PROC"

PDB="$1"; shift
[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [PSF/PDB] [TRJ]"; \
                  echo "mkvmd> [TRJ] could be multiple files."; \
                  echo "mkvmd> Note that XST files must exist in the same directory as the DCD files"; \
                  echo "mkvmd> By default, the MHC is '$SEL_MHC', the peptide is '$SEL_PEP'"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

cat > tcl <<EOF
package require pbctools
mol new $PDB waitfor all
EOF

for ii in $@; do
    prefix=${ii%.*}
    echo "mol addfile $ii waitfor all" >> tcl
    echo "pbc readxst ${prefix}.xst" >> tcl
done

cat >> tcl <<EOF
pbc wrap -all -compound residue -center com -centersel "${SEL_MHC}"
set sel_all [atomselect top all]
set sel_ref0 [atomselect top "${SEL_MHC}" frame 0]
set sel_ref [atomselect top "${SEL_MHC}"]
set sel_rmsd0 [atomselect top "${SEL_PEP}" frame 0]
set sel_rmsd [atomselect top "${SEL_PEP}"]
set num_frames [molinfo top get numframes]
set outfile [open "peprmsd" "w"]
for {set ii 0} {\$ii<\$num_frames} {incr ii} {
    \$sel_all frame \$ii
    \$sel_ref frame \$ii
    \$sel_rmsd frame \$ii
    \$sel_all move [measure fit \$sel_ref \$sel_ref0]
    set rmsd [measure rmsd \$sel_rmsd \$sel_rmsd0]
    puts \$outfile "\$ii \t \$rmsd"
}
#animate write dcd wrapped.dcd
quit
EOF

vmd -dispdev text -e tcl
rm tcl
