set selp [atomselect 0 all]
set selmem [atomselect 1 all]
set A [measure center $selmem]
# $selmem moveby [vecinvert $A]
set B [vecsub [measure center $selmem] [measure center $selp]]
# $selp moveby $B

package require Orient
namespace import Orient::orient

set sel [atomselect 0 "protein and name CA"]

draw delete all
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 0] {0 0 -1}]
$selp move $A
draw delete all
set I [draw principalaxes $sel]
# set A [orient $sel [lindex $I 0] {1 0 0}]
# $sel1 move $A
# draw delete all
# set I [draw principalaxes $sel]