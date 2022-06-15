#!/bin/bash
VMD="/opt/vmd/1.9.3/vmd -dispdev text"
NAMD="namd3 +p8 +devices 0"

[ $# -eq 0 ] && { echo "mknamd> Usage: $0 [-norun|-run]"; echo "mknamd> "; exit 1; }


# /-------------------/
# /     Functions     /
# /-------------------/

psfgen_header () {
  cat > tcl <<'EOF'
package require psfgen

mol new pdb2namd/md.pdb
foreach ii {A B} {
    set sel [atomselect top "segname PRO$ii"]
    $sel writepdb chains/PRO$ii.pdb
}
set sel [atomselect top "segname PROC and resid 1 2 3 4 5 6 7 8 9 and not name C CA N O HN HA CB"]
foreach name [$sel get name] {
  set sel [atomselect top "segname PROC and resid 1 2 3 4 5 6 7 8 9 and name $name"]
  $sel set name ${name}A
}
set sel [atomselect top "segname PROC"]
$sel writepdb chains/mutant.pdb

resetpsf
topology pdb2namd/readcharmmtop1.2/top_all36_prot.rtf
topology pdb2namd/readcharmmtop1.2/top_all36_hybrid.inp
# topology pdb2namd/top_all36_propatch.rtf
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
  mutate 1 2G
  mutate 2 2I
  mutate 3 2L
  mutate 4 2G
  mutate 5 2F
  mutate 6 2V
  mutate 7 2F
  mutate 8 2T
  mutate 9 2L
  first none
}
# patch AABP MUT:5
# patch AASP MUT:6
coordpdb chains/mutant.pdb MUT

EOF
}

psfgen_free () {
  cat >> tcl <<'EOF'
regenerate angles dihedrals
guesscoord
writepsf prot.psf
writepdb prot.pdb
quit
EOF
  $VMD -dispdev text -e tcl >& LOG_vmd
}

psfgen_bound () {
  cat >> tcl <<'EOF'
segment PROA {
  pdb chains/PROA.pdb
}
patch DISU PROA:101 PROA:164
patch DISU PROA:203 PROA:259
coordpdb chains/PROA.pdb PROA

segment PROB {
  pdb chains/PROB.pdb
}
patch DISU PROB:25 PROB:80
coordpdb chains/PROB.pdb PROB

regenerate angles dihedrals
guesscoord
writepsf prot.psf
writepdb prot.pdb
quit
EOF
  $VMD -dispdev text -e tcl >& LOG_vmd
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
solvate $filename.psf tmp.pdb -minmax {{-38 -40 -51} {38 40 51}} -o solvated
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
  $VMD -dispdev text -e tcl2 >& /dev/null
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
  $VMD -dispdev text -e tcl3 >& /dev/null
}


fixbb () {
  cat > tcl4 <<EOF
mol new pdb2namd/vmd_solvate/ionized.pdb type pdb waitfor all
set all [atomselect top "all"]
\$all set beta 0
set sel [atomselect top "segname MUT and name C N O"]
\$sel get resid
\$sel get resname
\$sel get name
\$sel set beta 1
\$all writepdb cons.fep
quit
EOF
  $VMD -dispdev text -e tcl4 >& /dev/null
}




# /-------------------/
# /     Main body     /
# /-------------------/

mkdir -p chains

echo "mknamd> Preparing Bound State"
[ -d bound ] && { cp -r bound bound.BAK; rm -r bound; }
mkdir -p bound bound/pdb2namd bound/pdb2namd/vmd_solvate

# /---------------------/
# /     Bound State     /
# /---------------------/

psfgen_header
psfgen_bound
solvate
markfep
fixbb

mv *.fep bound
mv ionized.p* cell_size.str bound/pdb2namd/vmd_solvate
cd bound
rsync -rpt ../eq/fep.tcl eq/
rsync -rpt ../eq/fep.eq.namd eq/fep.eq.namd
rsync -rpt ../eq/fep.namd eq/fep.namd
rsync -rpt ../mknamd_fep_check.sh mknamd_fep_check.sh
rsync -rpt ../mknamd_fep_run.sh mknamd_fep_run.sh
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

psfgen_header
psfgen_free
solvate
markfep
fixbb

mv *.fep free
mv ionized.p* cell_size.str free/pdb2namd/vmd_solvate
cd free
rsync -rpt ../eq/fep.tcl eq/
rsync -rpt ../eq/fep.eq.namd eq/fep.eq.namd
rsync -rpt ../eq/fep.namd eq/fep.namd
rsync -rpt ../mknamd_fep_check.sh mknamd_fep_check.sh
rsync -rpt ../mknamd_fep_run.sh mknamd_fep_run.sh
ln -s /home/kevin/forcefield/namd/charmm36/toppar_c36_jul20/ toppar
ln -s /home/kevin/forcefield/namd/charmm36/toppar_water_ions_namd.str .

if [[ "$1" == "-run" ]]; then 
  cd eq
  echo "mknamd> Running eq..."
  $NAMD fep.eq.namd >& LOG_eq
  cd ..
fi
cd ..

rm -r chains 
rm tcl* tmp.pdb prot.p* solvated.*
