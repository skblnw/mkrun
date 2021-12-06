#!/bin/bash
VMD="/opt/vmd/1.9.3/vmd"

SELECT_TEXT="segname RBD"

[ $# -eq 0 ] && { echo "mknamd> Usage: $0 [-norun|-run]"; echo "mknamd> Default peptide selection is: $SELECT_TEXT"; echo "mknamd> If apply multiple RESIDs, use e.g. \"1 2 3\""; exit 1; }

# length_of_peptide=`grep PROC pdb2namd/md.pdb | grep "CA" -c`
resid=(339 371 373 375 417 440 446 477 478 484 493 496 498 501 505)
sequence="GLY SER SER SER LYS ASN GLY SER THR GLU GLN GLY GLN ASN TYR"
# sequence=`grep PROC pdb2namd/md.pdb | grep "CA" | awk '{print $4}'`

declare -A mutation=(["ARG"]="R2A" \
                      ["ASN"]="N2A" \
                      ["ASP"]="D2A" \
                      ["CYS"]="C2A" \
                      ["GLN"]="Q2A" \
                      ["GLU"]="E2A" \
                      ["GLY"]="G2A" \
                      ["HIS"]="H2A" \
                      ["HSD"]="H2A" \
                      ["HSE"]="H2A" \
                      ["HSP"]="H2A" \
                      ["ILE"]="I2A" \
                      ["LEU"]="L2A" \
                      ["LYS"]="K2A" \
                      ["MET"]="M2A" \
                      ["PHE"]="F2A" \
                      ["PRO"]="P2A" \
                      ["SER"]="S2A" \
                      ["THR"]="T2A" \
                      ["TRP"]="W2A" \
                      ["TYR"]="Y2A" \
                      ["VAL"]="V2A" \
                      )


# /-------------------/
# /     Functions     /
# /-------------------/

psfgen () {
  cat > tcl <<EOF
package require psfgen 1.6

mol new pdb2namd/md.pdb
set sel [atomselect top "$SELECT_TEXT and resid $1 and not name C CA N O HN HA CB"]
foreach name [\$sel get name] {
  set sel [atomselect top "$SELECT_TEXT and resid $1 and name \$name"]
  \$sel set name \${name}A
}
set sel [atomselect top "$SELECT_TEXT"]
\$sel writepdb chains/mutant.pdb
EOF
  cat >> tcl <<'EOF'
resetpsf
topology pdb2namd/top_all36_propatch.rtf
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

EOF

# /--------------------------------/
# /     Check if it's Proline      /
# /     If yes, apply patches      /
# /--------------------------------/

  if [[ $resname == "PRO" ]]; then
    jj=$((ii-1))
    cat >> tcl <<EOF
segment MUT {
  pdb chains/mutant.pdb
  mutate $ii ${mutation[$resname]}
}
patch DISU MUT:336 MUT:361
patch DISU MUT:379 MUT:432
patch DISU MUT:391 MUT:525
patch DISU MUT:480 MUT:488
patch AABP MUT:$jj
patch AASP MUT:$ii
coordpdb chains/mutant.pdb MUT
EOF
  else
    cat >> tcl <<EOF
segment MUT {
  pdb chains/mutant.pdb
  mutate $ii ${mutation[$resname]}
}
patch DISU MUT:336 MUT:361
patch DISU MUT:379 MUT:432
patch DISU MUT:391 MUT:525
patch DISU MUT:480 MUT:488
coordpdb chains/mutant.pdb MUT
EOF
  fi

  cat >> tcl <<'EOF'
regenerate angles dihedrals
guesscoord
writepsf mutant.psf
writepdb mutant.pdb
quit
EOF
  $VMD -dispdev text -e tcl >& /dev/null
}

psfmerge () {
  cat > tcl2 <<EOF
package require topotools
mol new mutant.psf
mol addfile mutant.pdb
mol new pdb2namd/md.psf
mol addfile pdb2namd/md.pdb
set sel1 [atomselect 0 all]
set sel2 [atomselect 1 "segname PROA HETA HETB"]
set mol [::TopoTools::selections2mol "\$sel1 \$sel2"]
animate write psf prot.psf \$mol
animate write pdb prot.pdb \$mol
quit
EOF
  $VMD -dispdev text -e tcl2 >& /dev/null
}

solvate () {
  cat > tcl3 <<'EOF'
package require Orient
mol new prot.pdb
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
#solvate prot.psf tmp.pdb -t 10 -o solvated
solvate prot.psf tmp.pdb -minmax {{-43 -49 -70} {43 49 70}} -o solvated
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
  $VMD -dispdev text -e tcl3 >& /dev/null
}

markfep () {
  cat > tcl4 <<EOF
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
  $VMD -dispdev text -e tcl4 >& /dev/null
}







# /-------------------/
# /     Main body     /
# /-------------------/

nn=0
for resname in $sequence
do
  ii=${resid[$nn]}
  echo "mknamd> $ii $resname"

  [ -d posi$ii ] && { cp -r posi$ii posi$ii.BAK; rm -r posi$ii; }
  mkdir -p chains posi$ii 

# /---------------------/
# /     Bound State     /
# /---------------------/
  # [ -d posi$ii/bound ] && { cp -r posi$ii/bound posi$ii/bound.BAK; rm -r posi$ii/bound; }

  psfgen $ii
  psfmerge
  solvate
  markfep

  mkdir -p posi$ii/bound posi$ii/bound/pdb2namd posi$ii/bound/pdb2namd/vmd_solvate
  mv ionized.fep posi$ii/bound
  mv ionized.p* cell_size.str posi$ii/bound/pdb2namd/vmd_solvate
  cd posi$ii/bound
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.tcl eq/fep.tcl
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.eq.namd eq/fep.eq.namd
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.namd eq/fep.namd
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/mknamd_fep_check.sh mknamd_fep_check.sh
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/mknamd_fep_run.sh mknamd_fep_run.sh
  ln -s /home/kevin/forcefield/namd/charmm36/toppar_c36_jul20/ toppar
  ln -s /home/kevin/forcefield/namd/charmm36/toppar_water_ions_namd.str .

  if [[ "$1" == "-run" ]]; then 
    cd eq
    echo "mknamd> Running eq..."
    namd3 +p8 +devices 0 fep.eq.namd >& LOG_eq
    cd ..
  fi
  cd ../..

# /--------------------/
# /     Free State     /
# /--------------------/
  # [ -d posi$ii/free ] && { cp -r posi$ii/free posi$ii/free.BAK; rm -r posi$ii/free; }

  psfgen $ii
  mv mutant.psf prot.psf
  mv mutant.pdb prot.pdb
  solvate
  markfep

  mkdir -p posi$ii/free posi$ii/free/pdb2namd posi$ii/free/pdb2namd/vmd_solvate
  mv ionized.fep posi$ii/free
  mv ionized.p* cell_size.str posi$ii/free/pdb2namd/vmd_solvate
  cd posi$ii/free
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.tcl eq/fep.tcl
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.eq.namd eq/fep.eq.namd
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/eq/fep.namd eq/fep.namd
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/mknamd_fep_check.sh mknamd_fep_check.sh
  rsync -rpt /home/kevin/github/mkrun/NAMD/fep/mknamd_fep_run.sh mknamd_fep_run.sh
  ln -s /home/kevin/forcefield/namd/charmm36/toppar_c36_jul20/ toppar
  ln -s /home/kevin/forcefield/namd/charmm36/toppar_water_ions_namd.str .

  if [[ "$1" == "-run" ]]; then 
    cd eq
    echo "mknamd> Running eq..."
    namd3 +p8 +devices 0 fep.eq.namd >& LOG_eq
    cd ..
  fi
  cd ../..

  rm -rf chains tcl* tmp.pdb prot.p* solvated.* mutant.p*

  res_prev=$resname
  nn=$((nn+1))
done
