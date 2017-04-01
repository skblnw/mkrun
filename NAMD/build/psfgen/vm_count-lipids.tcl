set FILENAME membrane_psfgen
mol new $FILENAME.psf
mol addfile $FILENAME.pdb

set selall [atomselect top all]
set selmem [atomselect top "not {protein or water or ion}"]
#set selmemp [atomselect top "segname O1 and name P P1 P2 P3"]
set center [measure center $selmem]; list
set zcenter [lindex $center 2]; list

puts "In the membrane, I found:"
set reslist [lsort -unique [$selmem get resname]]

puts "In total, there are:"
set selup [atomselect top "not {protein or water or ion} and element P and z > ${zcenter}"]; list
puts "\[UP\]    [llength [lsort -unique [$selup get residue]]]"
set seldown [atomselect top "not {protein or water or ion} and element P and z < ${zcenter}"]; list
puts "\[DOWN\]  [llength [lsort -unique [$seldown get residue]]]"

foreach res $reslist {
    puts "${res}:"
    set selup [atomselect top "not {protein or water or ion} and resname ${res} and element P and z > ${zcenter}"]
    puts [llength [lsort -unique [$selup get residue]]]
    set seldown [atomselect top "not {protein or water or ion} and resname ${res} and element P and z < ${zcenter}"]
    puts [llength [lsort -unique [$seldown get residue]]]
    $selup delete
    $seldown delete
}

quit