#!/bin/bash
VMD="/opt/vmd/1.9.3/vmd"

rm -rf chains prot*
mkdir chains

cat > tcl <<'EOF'
package require psfgen 1.6

mol new md.pdb
foreach ii {A B C} {
    set sel [atomselect top "chain $ii"]
    $sel writepdb chains/$ii.pdb
}

resetpsf
topology readcharmmtop1.2/top_all36_prot.rtf
topology readcharmmtop1.2/top_all36_hybrid.inp
topology readcharmmtop1.2/top_all27_prot_lipid_na.inp
#topology readcharmmtop1.2/toppar_water_ions_namd.str

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

foreach ii {B C} {
  segment PRO$ii {pdb chains/$ii.pdb}
  coordpdb chains/$ii.pdb PRO$ii
}

foreach ii {A} {
  segment MUT {
    pdb chains/$ii.pdb
    mutate 478 T2K
  }
  patch DISU MUT:480 MUT:488
  coordpdb chains/$ii.pdb MUT
}

regenerate angles dihedrals
guesscoord
writepsf prot.psf
writepdb prot.pdb
quit
EOF
$VMD -dispdev text -e tcl 
rm -f tcl tmp.p*
