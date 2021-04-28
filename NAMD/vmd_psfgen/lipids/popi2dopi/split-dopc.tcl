mol new seg-DOPC.pdb

set dopi [atomselect top "residue 0 to 15 or residue 80 to 100"]
$dopi writepdb seg-DOPI.pdb

set dopi [atomselect top "not {residue 0 to 15 or residue 80 to 100}"]
$dopi writepdb seg-DOPC-new.pdb


mol delete all
mol new seg-DOPI.pdb
mol new seg-DOPC-new.pdb
