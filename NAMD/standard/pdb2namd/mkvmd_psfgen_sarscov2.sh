#!/bin/bash
VMD="/opt/vmd/1.9.3/vmd"

rm -rf chains mutant.p* prot* 
mkdir chains
DIR="/data/kevin/sarscov2/raw/6m0j-charmm-gui-3798750918"

cat > tcl <<EOF
mol new $DIR/step1_pdbreader.psf
mol addfile $DIR/step1_pdbreader.pdb
set sel [atomselect top "segname PROB"]
\$sel writepdb chains/RBD.pdb
quit
EOF
$VMD -dispdev text -e tcl 

cat > tcl <<'EOF'
package require psfgen
resetpsf
paratypeCharmm on 
mergeCrossterms no
topology top_all36_prot.rtf
# topology toppar_water_ions_namd.str

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

segment RBD {
  pdb chains/RBD.pdb
  mutate 339 ASP
  mutate 371 LEU
  mutate 373 PRO
  mutate 375 PHE
  mutate 417 ASN
  mutate 440 LYS
  mutate 446 SER
  mutate 477 ASN
  mutate 478 LYS
  mutate 484 ALA
  mutate 493 ARG
  mutate 496 SER
  mutate 498 ARG
  mutate 501 TYR
  mutate 505 HSD
  first NTER
  last CTER
}
patch DISU RBD:336 RBD:361
patch DISU RBD:379 RBD:432
patch DISU RBD:391 RBD:525
patch DISU RBD:480 RBD:488
coordpdb chains/RBD.pdb RBD

guesscoord
writepsf mutant.psf
writepdb mutant.pdb
quit
EOF
$VMD -dispdev text -e tcl 

cat > tcl <<EOF
package require topotools
mol new mutant.psf
mol addfile mutant.pdb
mol new $DIR/step1_pdbreader.psf
mol addfile $DIR/step1_pdbreader.pdb
set sel1 [atomselect 0 all]
set sel2 [atomselect 1 "segname PROA HETA HETB WATA"]
set mol [::TopoTools::selections2mol "\$sel1 \$sel2"]
animate write psf prot.psf \$mol
animate write pdb prot.pdb \$mol
quit
EOF
$VMD -dispdev text -e tcl 
rm tcl
