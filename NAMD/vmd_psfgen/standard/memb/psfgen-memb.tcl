package require psfgen
resetpsf

set HOME_DIR /home/kevin/
set FFDIR $HOME_DIR/force-field/charmm36

# CHARMM36 topologies
topology ${FFDIR}/toppar/top_all36_carb.rtf
topology ${FFDIR}/toppar/top_all36_cgenff.rtf
topology ${FFDIR}/toppar/top_all36_lipid.rtf
topology ${FFDIR}/toppar/top_all36_na.rtf
topology ${FFDIR}/toppar/top_all36_prot.rtf
topology ${FFDIR}/toppar/toppar_all36_lipid_cardiolipin.str
topology ${FFDIR}/toppar/toppar_all36_lipid_inositol.str
# topology of water is specially prepared for NAMD
topology ${FFDIR}/toppar_water_ions_namd.str
# User-defined topologies
#topology ${FFDIR}/toppar_dopi.str
#topology ${FFDIR}/toppar_tocl2.str

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
