#!/bin/zsh
source ~/.zshrc

REF_RTF='drawing_3D/drawing_3D.rtf'
LIGAND_NAME='LIG'
LIGAND_CHARGE=0

cat > hybrid.rtf <<'EOF'
* ---- 
* Built RTF for drawing_3D.mol2 
*    by user vzoete     Wed Sep  8 13:04:31 UTC 2021 
* ---- 
*
   22    0

EOF
grep "^MASS" $REF_RTF >> hybrid.rtf
cat >> hybrid.rtf <<EOF

AUTOGENERATE ANGLES DIHE
DEFA FIRS NONE LAST NONE

RESI $LIGAND_NAME  $LIGAND_CHARGE
GROUP
EOF
grep "^ATOM" drawing_3D/drawing_3D.rtf | awk '{$2 = ($2 != "C" && 
                                                     $2 != "CA" && 
                                                     $2 != "O" && 
                                                     $2 != "OXT" && 
                                                     $2 != "N" && 
                                                     $2 != "HA" && 
                                                     $2 != "H1" && 
                                                     $2 != "H2" && 
                                                     $2 != "H" \
                                                     ? $2"A" : $2)} 1' | column -t >> hybrid.rtf

cat >> hybrid.rtf <<'EOF'
GROUP                   
ATOM CBB  CT3    -0.27  
ATOM HB1B HA3     0.09  
ATOM HB2B HA3     0.09  
ATOM HB3B HA3     0.09  
EOF

list=`grep "^ATOM" drawing_3D/drawing_3D.rtf | awk '$2 !~ /^C$|^CA$|^O$|^OXT$|^N$|^HA$|^H1$|^H2$|^H$/ {print $2}'`
grep "^BOND" drawing_3D/drawing_3D.rtf > tmp
for ii in $list; do
    sed -i '' 's/'$ii'/'$ii'A/g' tmp
done
column -t tmp >> hybrid.rtf 
cat >> hybrid.rtf <<'EOF'
BOND CBB CA  
BOND CBB HB1B 
BOND CBB HB2B 
BOND CBB HB3B 
EOF
grep "^IMPH" drawing_3D/drawing_3D.rtf > tmp
for ii in $list; do
    sed -i '' 's/'$ii'/'$ii'A/g' tmp
done
column -t tmp >> hybrid.rtf 
grep "^IC" drawing_3D/drawing_3D.rtf > tmp
for ii in $list; do
    sed -i '' 's/'$ii'/'$ii'A/g' tmp
done
column -t tmp >> hybrid.rtf  
cat >> hybrid.rtf <<'EOF'
IC N    C    *CA   CBB   1.4592 114.4400  123.2300 111.0900  1.5461
IC C    CA   CBB   HB1B  1.5390 111.0900  177.2500 109.6000  1.1109
IC HB1B CA   *CBB  HB2B  1.1109 109.6000  119.1300 111.0500  1.1119
IC HB1B CA   *CBB  HB3B  1.1109 109.6000 -119.5800 111.6100  1.1114
EOF



cat > tcl <<'EOF'
package require psfgen
resetpsf
topology readcharmmtop1.2/top_all36_prot.rtf
topology readcharmmtop1.2/top_all36_hybrid.inp
topology toppar_water_ions_namd.str
topology hybrid.rtf
segment LIG {
    pdb drawing_3D/drawing_3D.pdb
    mutate 1 LIG}
coordpdb drawing_3D/drawing_3D.pdb LIG
regenerate angles dihedrals
guesscoord
writepsf prot.psf
writepdb prot.pdb
quit
EOF

vmd -dispdev text -e tcl
rm tcl tmp