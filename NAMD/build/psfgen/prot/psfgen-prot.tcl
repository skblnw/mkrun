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
# Topology of water is specially prepared for NAMD
topology ${FFDIR}/toppar_water_ions_namd.str
# User-defined topologies

# Aliases borrowed from AutoPSF
  pdbalias residue G GUA
  pdbalias residue C CYT
  pdbalias residue A ADE
  pdbalias residue T THY
  pdbalias residue U URA

  foreach bp { GUA CYT ADE THY URA } {
     pdbalias atom $bp "O5\*" O5'
     pdbalias atom $bp "C5\*" C5'
     pdbalias atom $bp "O4\*" O4'
     pdbalias atom $bp "C4\*" C4'
     pdbalias atom $bp "C3\*" C3'
     pdbalias atom $bp "O3\*" O3'
     pdbalias atom $bp "C2\*" C2'
     pdbalias atom $bp "O2\*" O2'
     pdbalias atom $bp "C1\*" C1'
  }

  pdbalias atom ILE CD1 CD
  pdbalias atom SER HG HG1
  pdbalias residue HIS HSD

# Heme aliases
  pdbalias residue HEM HEME
  pdbalias atom HEME "N A" NA
  pdbalias atom HEME "N B" NB
  pdbalias atom HEME "N C" NC
  pdbalias atom HEME "N D" ND

# Water aliases
  pdbalias residue HOH TIP3
  pdbalias atom TIP3 O OH2

# Ion aliases
  pdbalias residue K POT
  pdbalias atom K K POT
  pdbalias residue ICL CLA
  pdbalias atom ICL CL CLA
  pdbalias residue INA SOD
  pdbalias atom INA NA SOD
  pdbalias residue CA CAL
  pdbalias atom CA CA CAL

# Other aliases
  pdbalias atom LYS 1HZ HZ1
  pdbalias atom LYS 2HZ HZ2
  pdbalias atom LYS 3HZ HZ3

  pdbalias atom ARG 1HH1 HH11
  pdbalias atom ARG 2HH1 HH12
  pdbalias atom ARG 1HH2 HH21
  pdbalias atom ARG 2HH2 HH22

  pdbalias atom ASN 1HD2 HD21
  pdbalias atom ASN 2HD2 HD22


# User-defined alias
  pdbalias residue MSE MET

# Configure chains
foreach seg {P1 W1 P2 W2} {
    segment $seg {pdb seg-$seg.pdb}
    coordpdb seg-$seg.pdb $seg
}

guesscoord
writepdb prot_psfgen.pdb
writepsf prot_psfgen.psf

quit
