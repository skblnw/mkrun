mol new 1uru_dimer.pdb

set nn 0
foreach chain {A B} {
    incr nn 1
    foreach {atmtype segment} {protein P water W} {
        set sel [atomselect top "chain $chain and $atmtype"]
        $sel set segname ${segment}${nn}
        $sel set chain ${segment}
        if {0} {
            set oldNumbering [$sel get residue]
            set newNumbering {}
            foreach resid $oldNumbering {
                lappend newNumbering [expr {$resid + 1}]
            }
            $sel set resid $newNumbering
        }
    $sel writepdb seg-${segment}${nn}.pdb
    }    
}
