package require psfgen
resetpsf

set FFDIR /home/kevin/Dropbox/QWD/force-field/charmm36

topology ${FFDIR}/toppar/top_all36_lipid.rtf
topology ${FFDIR}/toppar/top_all36_carb.rtf
topology ${FFDIR}/toppar/toppar_all36_lipid_cardiolipin.str
topology ${FFDIR}/toppar/toppar_all36_lipid_inositol.str
topology ${FFDIR}/toppar_dopi.str
topology ${FFDIR}/toppar_tocl2.str

pdbalias residue POPI DOPI

foreach lipid {DOPI} {
  segment $lipid {pdb seg-POPI.pdb}
  coordpdb seg-POPI.pdb $lipid
}

foreach lipid {DOPC DOPE DOPA DOPS DOPG} {
  segment $lipid {pdb seg-$lipid.pdb}
  coordpdb seg-$lipid.pdb $lipid
}

segment TOCL {pdb seg-TOCL2.pdb}
coordpdb seg-TOCL2.pdb TOCL

guesscoord
writepdb membrane_psfgen.pdb
writepsf membrane_psfgen.psf

quit
