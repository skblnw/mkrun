mol new test.psf
mol addfile test.pdb

foreach {lipid resid} {DOPC 34 DOPI 2} {
set sel [atomselect top "resname $lipid and resid $resid"]
$sel writepdb 1-$lipid.pdb
}


quit
