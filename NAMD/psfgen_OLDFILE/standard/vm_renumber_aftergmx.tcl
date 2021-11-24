mol new seg_50mol.pdb

set selall [atomselect top all]
$selall set chain P

set NSTART 0
for {set ii 1} {$ii <= 50} {incr ii} {
    set NEND [expr $NSTART + 90]

    set sel [atomselect top "index $NSTART to $NEND"]
    $sel set segname P$ii

    set oldNumbering [$sel get residue]
    set newNumbering {}
    foreach resid $oldNumbering {
        lappend newNumbering [expr {$resid % 6 + 1}]
    }
    $sel set resid $newNumbering

    $sel writepdb P$ii.pdb
    set NSTART [expr $NEND + 1]
}

#$selall writepdb all.pdb

quit
