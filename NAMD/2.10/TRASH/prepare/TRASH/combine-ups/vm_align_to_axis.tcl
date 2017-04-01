package require Orient
namespace import Orient::orient

set sel [atomselect 0 all]
#set selopen [atomselect top "segname P1 and {resid 61 to 74 or resid 130 to 135 or resid 100 to 105}"]
set selopen [atomselect 0 "protein and name CA"]

draw delete all
set I [draw principalaxes $selopen]
set A [orient $selopen [lindex $I 0] {0 0 -1}]
$sel move $A
draw delete all
set I [draw principalaxes $selopen]