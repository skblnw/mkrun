proc get_cell {} { 
    set fout [open "cell_size.str" w]

    set all [atomselect top water] 
    set minmax [measure minmax $all] 
    set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]] 
    puts $fout "cellBasisVector1 [lindex $vec 0] 0 0" 
    puts $fout "cellBasisVector2 0 [lindex $vec 1] 0" 
    puts $fout "cellBasisVector3 0 0 [lindex $vec 2]" 
    set center [measure center $all] 
    puts $fout "cellOrigin $center" 

    $all delete 
    close $fout
}

mol new ionized.pdb
get_cell
quit
