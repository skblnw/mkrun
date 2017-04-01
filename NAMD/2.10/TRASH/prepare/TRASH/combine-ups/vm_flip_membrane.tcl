package require Orient
namespace import Orient::orient

set sel [atomselect top "name P"]
set selall [atomselect top "all"]


set I0 [draw principalaxes $sel]
set A [orient $sel [lindex $I0 0] {0 0 -1}]
$selall move $A
set I1 [draw principalaxes $sel]
# set A [orient $sel [lindex $I1 1] [vecinvert [lindex $I0 1]]]
# $selall move $A
# set I1 [draw principalaxes $sel]