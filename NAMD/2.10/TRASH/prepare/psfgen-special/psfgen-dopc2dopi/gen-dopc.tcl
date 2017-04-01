package require psfgen

topology toppar/top_all36_lipid.rtf
topology toppar/top_all36_carb.rtf
topology toppar/toppar_all36_lipid_cardiolipin.str
topology toppar/toppar_all36_lipid_inositol.str

segment DOPC {pdb seg-DOPC.pdb}

coordpdb seg-DOPC.pdb DOPC

guesscoord
writepdb seg-DOPC-old.pdb
writepsf seg-DOPC-old.psf

quit
