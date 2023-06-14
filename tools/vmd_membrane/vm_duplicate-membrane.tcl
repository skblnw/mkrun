#########################################
## Description: Duplicate small membranes to create large ones
## Author: Kev Dec 2016
#########################################

# Assuming your membrane is square in shape
# Usually it requires you to run an equilibrated run of the membrane first
# Then from the xsc file, you will have the boxsize (same for x and y)
set boxsize 51.0363121466
# Set the final size of the duplication
# It will affect the times of iterations along each side
set xx 4
set yy 3

set filename "../charmm-gui/step4_lipid"

set kk 1
set arg_sel ""
for {set ii 0} {$ii < $xx} {incr ii} {
    for {set jj 0} {$jj < $yy} {incr jj} {
        mol new $filename.psf waitfor all
        mol addfile $filename.pdb waitfor all
        set sel_index [set sel($kk) [atomselect top all]]
        lappend arg_sel "$sel_index"
        $sel($kk) set segname [join "L$kk"]

        set tmp [expr $boxsize * $ii]
        set dist [list $tmp 0 0]
        $sel($kk) moveby $dist
        
        set tmp [expr $boxsize * $jj]
        set dist [list 0 $tmp 0]
        $sel($kk) moveby $dist

        set kk [incr kk]
    }
}

package require topotools
puts $arg_sel
set mol [::TopoTools::selections2mol "$arg_sel"]

animate write psf duplicated.psf $mol
animate write pdb duplicated.pdb $mol
quit
