mol new membrane-charmmgui.pdb

foreach lipid {DOPC DOPE DOPA DOPG DOPS TOCL} {
    set seg [atomselect top "resname $lipid"]
    $seg writepdb seg-$lipid.pdb
}

mol delete all
foreach lipid {DOPC DOPE DOPA DOPG DOPS TOCL} {
    mol new seg-$lipid.pdb
}