
set prefix "t1.400ns"

proc writeXSC { prefix } {
  set outXSC [open "${prefix}.restart.xsc" w]
  puts $outXSC "\# NAMD extended system configuration output file"
  puts $outXSC "\#\\\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
  puts $outXSC "0 [molinfo top get a] 0 0 0 [molinfo top get b] 0 0 0 [molinfo top get c] [lindex [molinfo top get center] 0 0] [lindex [molinfo top get center] 0 1] [lindex [molinfo top get center] 0 2] 0 0 0 0 0 0"
  close $outXSC
}

set all [atomselect top all]
$all frame now
$all writenamdbin $prefix.restart.coor
writeXSC $prefix
