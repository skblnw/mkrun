#!/bin/bash

cat > tcl <<EOF
set outf [open tmdlist.txt "w"]
EOF

for ii in $(seq 1 5); do
    cat >> tcl <<EOF
mol new ${ii}a.pdb
set sel [atomselect top "beta > 0"]
set outline [\$sel get resid]
mol new ${ii}b.pdb
set sel [atomselect top "beta > 0"]
foreach xx [\$sel get resid] {lappend outline \$xx}
puts \$outf \$outline
EOF
done

cat >> tcl <<EOF
close \$outf
quit
EOF

vmd -dispdev text -e tcl > /dev/null


cat > tcl <<EOF
set outf [open reslist.txt "w"]
EOF

for ii in $(seq 1 14); do
    cat >> tcl <<EOF
mol new res${ii}a.pdb
set sel [atomselect top "beta > 0"]
set outline [\$sel get resid]
mol new res${ii}b.pdb
set sel [atomselect top "beta > 0"]
foreach xx [\$sel get resid] {lappend outline \$xx}
puts \$outf \$outline
EOF
done

cat >> tcl <<EOF
close \$outf
quit
EOF

vmd -dispdev text -e tcl > /dev/null
rm tcl

