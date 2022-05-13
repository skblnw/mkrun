#!/bin/bash

REF_RTF='swiss/drawing_3D.rtf'
LIGAND_NAME='LIG'
LIGAND_CHARGE=0
LIST="C4 C8 C9 H8"

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
grep "^ATOM" $REF_RTF | awk '{$2 = ($2 != "C4" && 
                                    $2 != "C8" &&
                                    $2 != "C9" &&  
                                    $2 != "H8" \
                                    ? $2 : $2"A")} 1' | column -t > tmp
cat >> tmp <<'EOF'
ATOM  C4B   CB    -0.1500
ATOM  H8B   HCMM   0.1500
EOF

column -t tmp >> hybrid.rtf 
list=`grep "^ATOM" $REF_RTF | awk '$2 !~ /^C4$|^C8$|^C9$|^H8$/ {print $2}'`
#################################################

grep "^BOND" $REF_RTF > tmp
for ii in $LIST; do
    sed -i 's/ '$ii' / '$ii'A/g' tmp
done
cat >> tmp <<'EOF'
BOND  C4B  H8B
BOND  C3   C4B  
BOND  C4B   C5  
EOF

column -t tmp >> hybrid.rtf 
#################################################

grep "^IMPH" $REF_RTF > tmp
for ii in $LIST; do
    sed -i 's/ '$ii' / '$ii'A/g' tmp
done
cat >> tmp <<'EOF'
IMPH C3   C4B   C2   H5   
IMPH C4B   C5   C3   C8   
IMPH C5   O1   C4B   C6
IMPH C4B   C3   C5   H8B
EOF

column -t tmp >> hybrid.rtf 
#################################################

grep "^IC" $REF_RTF > tmp
for ii in $LIST; do
    sed -i 's/ '$ii' / '$ii'A/g' tmp
done
cat >> tmp <<'EOF'
IC C2   C3    C4B   C5    1.39   120.00    0.01  120.03  1.39  
IC C2   C3    C4B   C8    1.39   120.00  179.99  120.03  1.29  
IC C3   C4B    C5   C6    1.39   120.03   -0.04  119.98  1.39  
IC C3   C4B    C5   O1    1.39   120.03 -179.95  120.00  1.35  
IC C3   C4B    C8   C9    1.39   120.03   13.53  179.93  1.19  
IC C4B   C3    C2   C7    1.39   120.00    0.01  120.00  1.39  
IC C4B   C3    C2   C1    1.39   120.00  179.98  120.01  1.46  
IC C4B   C5    C6   C7    1.39   119.98    0.05  120.02  1.39  
IC C4B   C5    C6   H6    1.39   119.98 -179.99  120.03  1.02  
IC C4B   C5    O1   S     1.39   120.00  -90.06  109.42  1.69  
IC C4B   C8    C9   H8    1.29   179.93 -176.95  179.92  0.92  
IC C5   C4B    C3   H5    1.39   120.03  179.97  120.05  1.02  
IC C5   C4B    C8   C9    1.39   119.94 -166.49  179.93  1.19  
IC C6   C5    C4B   C8    1.39   119.98  179.98  119.94  1.29  
IC O1   C5    C4B   C8    1.35   120.00    0.07  119.94  1.29  
IC C8   C4B    C3   H5    1.29   120.03   -0.05  120.05  1.02  
IC C4B   C2   *C3   H5    0.00     0.00  180.00    0.00  0.00  
IC C5   C3   *C4B   C8    0.00     0.00  180.00    0.00  0.00  
IC O1   C4B   *C5   C6    0.00     0.00  180.00    0.00  0.00
IC  C2  C3   C4B  H8B   1.39   119.97  169.96  120.00  1.02  
IC  C6  C5   C4B  H8B   1.39   120.00 -179.96  119.97  1.02  
IC  O1  C5   C4B  H8B   1.35   120.04    0.01  119.97  1.02  
IC  H8B C4B   C3  H5    1.02   120.00   -0.09  120.01  1.02  
IC  C3  C5  *C4B  H8B   0.00     0.00  180.00    0.00  0.00  
EOF

column -t tmp >> hybrid.rtf 
#################################################
