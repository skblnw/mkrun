mol new membrane_psfgen.psf
mol addfile membrane_psfgen.pdb

set selall [atomselect top all]
set selmem [atomselect top "segname L1 L2"]
set selmemp [atomselect top "segname L1 L2 and name P P1 P2 P3"]
set center [measure center $selmem]; list
set zcenter [lindex $center 2]; list

puts "In total:"
set selup [atomselect top "segname L1 L2 and name P P1 P2 P3 and z > ${zcenter}"]; list
puts "\[UP\]   [$selup num]"
set seldown [atomselect top "segname L1 L2 and name P P1 P2 P3 and z < ${zcenter}"]; list
puts "\[DOWN\] [$seldown num]"

set reslist [$selmemp get resname]; list
foreach res {DOPA DOPC DOPE DOPG DOPS DOPI TOCL} {
    if {[lsearch $reslist $res] >= 0} {
        puts "${res} is found:"
        set selup [atomselect top "segname L1 L2 and resname ${res} and name P P1 P2 P3 and z > ${zcenter}"]
        puts [$selup num]
        set seldown [atomselect top "segname L1 L2 and resname ${res} and name P P1 P2 P3 and z < ${zcenter}"]
        puts [$selup num]
        $selup delete
        $seldown delete
    } else {
        puts "88"
    }
}
