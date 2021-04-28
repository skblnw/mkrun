mol new seg-DOPC-old.psf
mol addfile seg-DOPC-old.pdb

set dopi [atomselect top "residue 0 to 15 or residue 80 to 95"]
$dopi writepdb seg-DOPI.pdb

set dopi [atomselect top "not {residue 0 to 15 or residue 80 to 95}"]
$dopi writepdb seg-DOPC-new.pdb


mol delete all
mol new seg-DOPI.pdb
mol new seg-DOPC-new.pdb
