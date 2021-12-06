#!/bin/bash
VMD="/opt/vmd/1.9.3/vmd -dispdev text"
NAMD="namd3 +p8 +devices 0"

[ $# -eq 0 ] && { echo "mknamd> Usage: $0 [-norun|-run]"; echo "mknamd> "; exit 1; }


# /-------------------/
# /     Functions     /
# /-------------------/

psfgen () {
  cat > tcl <<'EOF'
package require psfgen 1.6

mol new pdb2namd/md.pdb
foreach ii {A} {
    set sel [atomselect top "chain $ii"]
    $sel writepdb chains/$ii.pdb
}
set sel [atomselect top "name ZN"]
$sel writepdb chains/Z.pdb
set sel [atomselect top "chain B"]
$sel writepdb chains/mutant.pdb

resetpsf
topology pdb2namd/top_all36_propatch.rtf
topology pdb2namd/top_all36_hybridgly.rtf
topology pdb2namd/readcharmmtop1.2/top_all36_prot.rtf
topology pdb2namd/readcharmmtop1.2/top_all36_hybrid.inp
topology pdb2namd/toppar_water_ions_namd.str

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

segment MUT {
  pdb chains/mutant.pdb
  mutate 339 G2D
  mutate 371 S2L
  mutate 373 S2P
  mutate 375 S2F
  mutate 417 K2N
  mutate 440 N2K
  mutate 446 G2S
  mutate 477 S2N
  mutate 478 T2K
  mutate 484 E2A
  mutate 493 Q2R
  mutate 496 G2S
  mutate 498 Q2R
  mutate 501 N2Y
  mutate 505 Y2H
  first NTER
  last CTER
}
patch DISU MUT:336 MUT:361
patch DISU MUT:379 MUT:432
patch DISU MUT:391 MUT:525
patch DISU MUT:480 MUT:488
patch AABP MUT:372
patch AASP MUT:373
coordpdb chains/mutant.pdb MUT

foreach ii {A Z} {segment PRO$ii {pdb chains/$ii.pdb}; coordpdb chains/$ii.pdb PRO$ii}

regenerate angles dihedrals
guesscoord
writepsf prot.psf
writepdb prot.pdb
quit
EOF
  $VMD -e tcl >& LOG_vmd
}

solvate () {
  cat > tcl2 <<'EOF'
set filename prot

package require Orient
mol new $filename.pdb
set sel [atomselect top "all"]
set I [Orient::calc_principalaxes $sel]
set A [Orient::orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [Orient::calc_principalaxes $sel]
set A [Orient::orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
$sel moveby [vecinvert [measure center $sel]]
$sel writepdb tmp.pdb

package require solvate
#solvate $filename.psf tmp.pdb -t 10 -o solvated
solvate $filename.psf tmp.pdb -minmax {{-47 -47 -70} {47 47 70}} -o solvated
package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -o ionized

mol new ionized.pdb type pdb waitfor all
set all [atomselect top "all"]
set fout [open "cell_size.str" w]
set all [atomselect top water] 
set minmax [measure minmax $all] 
set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]] 
puts $fout "cellBasisVector1 [lindex $vec 0] 0 0" 
puts $fout "cellBasisVector2 0 [lindex $vec 1] 0" 
puts $fout "cellBasisVector3 0 0 [lindex $vec 2]" 
set center [measure center $all] 
puts $fout "cellOrigin $center" 
close $fout
quit
EOF
  $VMD -e tcl2 >& /dev/null
}

markfep () {
  cat > tcl3 <<EOF
mol new ionized.pdb
set sel [atomselect top all]
\$sel set beta 0
set sel [atomselect top "resname \".*2.*\" and not name CA CB HA HB and name \".*A\""]
\$sel set beta -1
set sel [atomselect top "resname \".*2.*\" and not name CA CB HA HB and name \".*B\""]
\$sel set beta 1
set sel [atomselect top all]
\$sel writepdb ionized.fep
quit
EOF
  $VMD -e tcl3 >& /dev/null
}







# /-------------------/
# /     Main body     /
# /-------------------/


echo "mknamd> Preparing Bound State"
[ -d bound ] && { cp -r bound bound.BAK; rm -r bound; }
mkdir -p chains bound bound/pdb2namd bound/pdb2namd/vmd_solvate

# /---------------------/
# /     Bound State     /
# /---------------------/

psfgen
cp tcl tcl.bak
solvate
markfep

mv ionized.fep bound
mv ionized.p* cell_size.str bound/pdb2namd/vmd_solvate
cd bound
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.tcl eq/
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.eq.namd eq/fep.eq.namd
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.namd eq/fep.namd
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/mknamd_fep_check.sh mknamd_fep_check.sh
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/mknamd_fep_run.sh mknamd_fep_run.sh
ln -s /home/kevin/forcefield/namd/charmm36/toppar_c36_jul20/ toppar
ln -s /home/kevin/forcefield/namd/charmm36/toppar_water_ions_namd.str .

if [[ "$1" == "-run" ]]; then 
  cd eq
  echo "mknamd> Running eq..."
  $NAMD fep.eq.namd >& LOG_eq
  cd ..
fi
cd ..

# /--------------------/
# /     Free State     /
# /--------------------/

echo "mknamd> Preparing Free State"
[ -d free ] && { cp -r free free.BAK; rm -r free; }
mkdir -p free free/pdb2namd free/pdb2namd/vmd_solvate

sed -i "s/.*} {segment PRO.*//g" tcl
$VMD -e tcl >& /dev/null
solvate
markfep

mv ionized.fep free
mv ionized.p* cell_size.str free/pdb2namd/vmd_solvate
cd free
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.tcl eq/
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.eq.namd eq/fep.eq.namd
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.namd eq/fep.namd
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/mknamd_fep_check.sh mknamd_fep_check.sh
rsync -rpt /home/kevin/github/mkrun/NAMD/fep/mknamd_fep_run.sh mknamd_fep_run.sh
ln -s /home/kevin/forcefield/namd/charmm36/toppar_c36_jul20/ toppar
ln -s /home/kevin/forcefield/namd/charmm36/toppar_water_ions_namd.str .

if [[ "$1" == "-run" ]]; then 
  cd eq
  echo "mknamd> Running eq..."
  $NAMD fep.eq.namd >& LOG_eq
  cd ..
fi
cd ..

rm -r chains tcl* tmp.pdb prot.p* solvated.*
