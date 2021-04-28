package require psfgen
resetpsf

topology toppar/top_all36_lipid.rtf
topology toppar/top_all36_carb.rtf
topology toppar/toppar_all36_lipid_cardiolipin.str
topology toppar/toppar_all36_lipid_inositol.str
topology toppar_dopi.str
topology toppar_tocl2.str

#pdbalias residue TOCL TOCL2

foreach lipid {DOPC DOPI} {
  segment $lipid {pdb seg-$lipid-new.pdb}
  coordpdb seg-$lipid-new.pdb $lipid
}

foreach lipid {DOPE DOPA DOPS DOPG} {
  segment $lipid {pdb seg-$lipid.pdb}
  coordpdb seg-$lipid.pdb $lipid
}

segment TOCL {pdb seg-TOCL.pdb}
coordpdb seg-TOCL.pdb TOCL

guesscoord
writepdb membrane_psfgen.pdb
writepsf membrane_psfgen.psf

quit
