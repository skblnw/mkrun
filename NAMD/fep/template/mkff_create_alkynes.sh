#!/bin/bash

REF_RTF='swiss/drawing_3D.rtf'
LIGAND_NAME='LIG'
LIGAND_CHARGE=0
LIST="C6 H7"

cat > hybrid.rtf <<EOF
* ---- 
* Built RTF for drawing_3D.mol2 
*    by user vzoete     Wed Sep  8 13:04:31 UTC 2021 
* ---- 
*
   22    0

MASS   201 CB    12.011000
MASS   202 C=O   12.011000
MASS   203 CR    12.011000
MASS   204 OR    15.999400
MASS   205 O=C   15.999400
MASS   206 NR    14.006700
MASS   207 SO2   32.066000
MASS   208 O2CM  15.999400
MASS   209 F     18.998400
MASS   210 CSP   12.011000
MASS   211 HCMM   1.007940
MASS   212 HOCO   1.007940
MASS   213 HNR    1.007940

AUTOGENERATE ANGLES DIHE
DEFA FIRS NONE LAST NONE

RESI $LIGAND_NAME  $LIGAND_CHARGE
GROUP
EOF
grep "^ATOM" $REF_RTF | awk '{$2 = ($2 != "C6" && 
                                    $2 != "H7" \
                                    ? $2 : $2"A")} 1' | column -t > tmp
cat >> tmp <<'EOF'
GROUP
ATOM  C6B   CB     0.0730
ATOM  C8B   CSP   -0.0730
ATOM  C9B   CSP   -0.1770
ATOM  H7B   HCMM   0.1770
EOF

column -t tmp >> hybrid.rtf 
# list=`grep "^ATOM" $REF_RTF | awk '$2 !~ /^C4$|^C8$|^C9$|^H8$/ {print $2}'`
#################################################

grep "^BOND" $REF_RTF > tmp
for ii in $LIST; do
    sed -i 's/ '$ii' / '$ii'A/g' tmp
done
cat >> tmp <<'EOF'
BOND C5    C6B
BOND C6B   C7
BOND C6B   C8B  
BOND C8B   C9B  
BOND C9B   H7B  
EOF

column -t tmp >> hybrid.rtf 
#################################################

grep "^IMPH" $REF_RTF > tmp
for ii in $LIST; do
    sed -i 's/ '$ii' / '$ii'A/g' tmp
done
cat >> tmp <<'EOF'
IMPH C5   O1   C4   C6B   
IMPH C7   C6B   C2   H8
IMPH C6B   C7   C5   C8B
EOF

column -t tmp >> hybrid.rtf 
#################################################

grep "^IC" $REF_RTF > tmp
for ii in $LIST; do
    sed -i 's/ '$ii' / '$ii'A/g' tmp
done
cat >> tmp <<'EOF'
IC C2   C7    C6B   C5    1.39   120.00    0.01  120.03  1.39  
IC C2   C7    C6B   C8B    1.39   120.00  179.99  120.03  1.29  
IC C7   C6B    C5   C4    1.39   120.03   -0.04  119.98  1.39  
IC C7   C6B    C5   O1    1.39   120.03 -179.95  120.00  1.35  
IC C6B   C7    C2   C3    1.39   120.00    0.01  120.00  1.39  
IC C6B   C7    C2   C1    1.39   120.00  179.98  120.01  1.46  
IC C6B   C5    C4   C3    1.39   119.98    0.05  120.02  1.39  
IC C6B   C5    C4   H6    1.39   119.98 -179.99  120.03  1.02  
IC C6B   C5    O1   S     1.39   120.00  -90.06  109.42  1.69  
IC C5   C6B    C7   H8    1.39   120.03  179.97  120.05  1.02  
IC C4   C5    C6B   C8B    1.39   119.98  179.98  119.94  1.29  
IC O1   C5    C6B   C8B    1.35   120.00    0.07  119.94  1.29  
IC C8B   C6B    C7   H8    1.29   120.03   -0.05  120.05  1.02  
IC C6B   C2   *C7   H8    0.00     0.00  180.00    0.00  0.00  
IC C5   C7   *C6B   C8B    0.00     0.00  180.00    0.00  0.00  
IC O1   C6B   *C5   C4    0.00     0.00  180.00    0.00  0.00
IC C7   C6B    C8B   C9B    1.39   120.03   13.53  179.93  1.19  
IC C6B   C8B    C9B   H7B    1.29   179.93 -176.95  179.92  0.92  
IC C5   C6B    C8B   C9B    1.39   119.94 -166.49  179.93  1.19
EOF

column -t tmp >> hybrid.rtf 
#################################################
