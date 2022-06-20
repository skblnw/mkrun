#!/bin/bash

SEL_TEXT="segname RBD"

PDB="$1"; shift
[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [PSF/PDB] [TRJ]"; \
                  echo "mkvmd> [TRJ] could be multiple files."; \
                  echo "mkvmd> Note that XST files must exist in the same directory as the DCD files"; \
                  echo "mkvmd> By default, the selection is '$SEL_TEXT'"; exit 1; }

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
pbc wrap -all -compound segid -center com -centersel "${SEL_TEXT}"
set sel_all [atomselect top all]
set sel_ref0 [atomselect top "${SEL_TEXT}" frame 0]
set sel_ref [atomselect top "${SEL_TEXT}"]
set num_frames [molinfo top get numframes]
for {set ii 0} {\$ii<\$num_frames} {incr ii} {
    \$sel_all frame \$ii
    \$sel_ref frame \$ii
    \$sel_all move [measure fit \$sel_ref \$sel_ref0]
}
animate write dcd wrapped.dcd
quit
EOF

vmd -dispdev text -e tcl
rm tcl
