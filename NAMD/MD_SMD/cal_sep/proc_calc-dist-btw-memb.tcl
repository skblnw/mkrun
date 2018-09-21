# Calc function as it is named
proc calc_dist_btw_memb_and_res {nn} {
    # select residue CAs according to resid
    set sel [atomselect top "protein" frame $nn]
    
    # select z of com of membrane phosphates
    set selmem [atomselect top "not {protein or water or ion}"]
    set center [measure center $selmem weight mass]
    set zcenter [lindex $center 2]
    set selup [atomselect top "not {protein or water or ion} and element P and z < ${zcenter}"]

    set dist [lindex [vecsub [measure center $sel weight mass] [measure center $selup weight mass]] 2]
    set res [format "%.2f" $dist]

    $sel delete

    return $res
}
