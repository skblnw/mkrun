#/bin/bash

for ii in $(seq 1 32)
do
    mkdir $ii
    cd $ii

    cat > leapin << EOF
    source leaprc.protein.ff14SB
    source leaprc.RNA.OL3
    source leaprc.DNA.bsc1
    source leaprc.water.tip3p

    mol = loadpdb ../../pdb/crd.pdb.$ii
    solvatebox mol TIP3PBOX 10
    addions mol K+ 176 Cl- 80
    saveamberparm mol ionized.parm7 ionized.rst7
    quit
EOF

    tleap -f leapin
    rm -f leapin
    cd ..

done
