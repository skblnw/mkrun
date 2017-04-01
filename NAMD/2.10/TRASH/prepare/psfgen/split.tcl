mol new ../raw/step4_lipid.psf
mol addfile ../raw/step4_lipid.pdb

foreach lipid {DOPC DOPE DOPA DOPG DOPS TOCL2 POPI} {
    set seg [atomselect top "resname $lipid"]
    $seg writepdb seg-$lipid.pdb
}

mol delete all
foreach lipid {DOPC DOPE DOPA DOPG DOPS TOCL2 POPI} {
    mol new seg-$lipid.pdb
}
