proc writeXSC { j fr } {
  set outXSC [open "win${j}.restart.xsc" w]
  animate goto $fr
  puts $outXSC "\# NAMD extended system configuration output file"
  puts $outXSC "\#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
  puts $outXSC "0 [molinfo top get a] 0 0 0 [molinfo top get b] 0 0 0 [molinfo top get c] [lindex [molinfo top get center] 0 0] [lindex [molinfo top get center] 0 1] [lindex [molinfo top get center] 0 2] 0 0 0 0 0 0"
  close $outXSC
}

set skip 5

## load files here

mol new build/ionized.psf
mol addfile SMD1/AmtB-smd01A.dcd step $skip waitfor all
mol addfile SMD1/AmtB-smd01B.dcd step $skip waitfor all
mol addfile SMD1/AmtB-smd01C.dcd step $skip waitfor all
mol addfile SMD2/AmtB-smd02A.dcd step $skip waitfor all


set sel1 [atomselect top "serial 28892"]
set sel2 [atomselect top "name CA"]

set all [atomselect top all]

for {set i 0} {$i < [molinfo top get numframes]} {incr i 1} {

   $sel1 frame $i
   $sel2 frame $i
   $all frame $i

   set c1 [measure center $sel1]
   set c2 [measure center $sel2]

   set distZ [lindex [vecsub $c1 $c2] 2]

   for {set j -13} {$j <= 5} {incr j} {
## note that this tolerance may need to be changed to capture all windows
      if {[expr abs($distZ - $j)] < 0.25} {
         $all writenamdbin win${j}.restart.coor
	 writeXSC $j $i

      }
   }

}

quit

