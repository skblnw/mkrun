#!/bin/bash

SEL1="segname PROA PROB"
SEL2="segname PROC PROD"
START=36
END=48

[ $# -ne 2 ] && { echo "mknamd> Usage: $0 [PSF/PDB] [TRJ]"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

PDB="$1"
TRJ="$2"

cat > tcl <<'EOF'
proc writeXSC { jj fr } {
  set outXSC [open "win${jj}.restart.xsc" w]
  animate goto $fr
  puts $outXSC "\# NAMD extended system configuration output file"
  puts $outXSC "\#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
  puts $outXSC "0 [molinfo top get a] 0 0 0 [molinfo top get b] 0 0 0 [molinfo top get c] [lindex [molinfo top get center] 0 0] [lindex [molinfo top get center] 0 1] [lindex [molinfo top get center] 0 2] 0 0 0 0 0 0"
  close $outXSC
}
EOF

cat >> tcl <<EOF
## load files here
mol new pdb2namd/vmd_solvate/ionized.psf
mol addfile smd-4.dcd waitfor all

set sel1 [atomselect top "$SEL1"]
set sel2 [atomselect top "$SEL2"]
set all [atomselect top all]

for {set ii 0} {\$ii < [molinfo top get numframes]} {incr ii 1} {

   \$sel1 frame \$ii
   \$sel2 frame \$ii
   \$all frame \$ii

   # set c1 [measure center \$sel1]
   # set c2 [measure center \$sel2]
   # set distZ [lindex [vecsub \$c2 \$c1] 2]
   set distZ [vecdist [measure center \$sel2 weight mass] [measure center \$sel1 weight mass]]

   for {set jj $START} {\$jj <= $END} {incr jj} {
      ## note that this tolerance may need to be changed to capture all windows
      if {[expr abs(\$distZ - \$jj)] < 1} {
         \$all writenamdbin win\${jj}.restart.coor
	      writeXSC \$jj \$ii
         continue
      }
   }
}
quit
EOF

vmd -dispdev text -e tcl
rm tcl
