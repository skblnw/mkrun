mol new ionized.pdb

set allatoms [atomselect top all]
$allatoms set beta 0
$allatoms set occupancy 0

# select z of com of membrane phosphates
set selmem [atomselect top "not {protein or water or ion}"]
set center [measure center $selmem weight mass]
set zcenter [lindex $center 2]
set sellow [atomselect top "not {protein or water or ion} and name P and z < ${zcenter}"]
$sellow set beta 1
puts "Membrane COM: [measure center $sellow weight mass]"

$allatoms writepdb ref.pdb


$allatoms set beta 0
$allatoms set occupancy 0

set selpro [atomselect top "protein"]
puts "Protein COM: [measure center $selpro weight mass]"
$selpro set beta 1

$allatoms writepdb main.pdb
quit
