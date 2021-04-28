mol new ../raw/step4_lipid.psf
mol addfile ../raw/step4_lipid.pdb

foreach lipid {DOPC DOPE DOPA DOPG DOPS TOCL2 POPI} {
    set seg [atomselect top "resname $lipid"]
    $sel1 set segname L1
    $sel1 set chain L
    set oldNumbering [$sel1 get residue]
    set newNumbering {}
    foreach resid $oldNumbering {
        lappend newNumbering [expr {$resid + 1}]
    }
    $sel1 set resid $newNumbering
    
    $seg writepdb seg-$lipid.pdb
}

mol delete all
foreach lipid {DOPC DOPE DOPA DOPG DOPS TOCL2 POPI} {
    mol new seg-$lipid.pdb
}