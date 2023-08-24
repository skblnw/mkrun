foreach ii {B C D} jj {B C D} {
    set sel [atomselect top "chain $ii and not {water or ion}"]
    $sel set chain $jj
    $sel writepdb $ii.pdb
}

foreach ii {MG NA HOH} jj {E F G} {
    set sel [atomselect top "resname $ii"]
    $sel set chain $jj
    $sel writepdb $jj.pdb
}