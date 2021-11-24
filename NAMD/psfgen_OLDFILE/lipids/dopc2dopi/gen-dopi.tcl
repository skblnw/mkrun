package require psfgen

topology toppar/top_all36_lipid.rtf
topology toppar/top_all36_carb.rtf
topology toppar/toppar_all36_lipid_cardiolipin.str
topology toppar/toppar_all36_lipid_inositol.str
topology toppar_dopi.str

pdbalias residue DOPC DOPI

segment DOPI {pdb seg-DOPI.pdb}

coordpdb seg-DOPI.pdb DOPI

guesscoord
writepdb seg-DOPI-new.pdb
writepsf seg-DOPI-new.psf

quit
