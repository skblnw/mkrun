mol new ionized.pdb

set allatoms [atomselect top all]
$allatoms set beta 0
$allatoms set occupancy 0

set selmem [atomselect top "not {water or ion} and not protein and name P and z < -30"]
set mempos [measure center $selmem]
$selmem set beta 1

$allatoms writepdb ref.pdb

$allatoms set beta 0
$allatoms set occupancy 0

set selpro [atomselect top "not {water or ion} and protein and serial 9941 to 11725 and name CA"]
set propos [measure center $selpro]
$selpro set beta 1

$allatoms writepdb main.pdb
quit
