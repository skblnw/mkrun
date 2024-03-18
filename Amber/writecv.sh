#!/bin/bash

PDB=$1

[ $# -ne 1 ] && { echo -e "Usage: $0 [PDB]"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

cat > tcl << EOF
mol new $PDB
EOF
cat >> tcl << 'EOF'
set outfile [open "cv.in" w]
set sel [atomselect top "backbone and resid 345 to 370"]
puts $outfile " &colvar"
puts $outfile "    cv_type = 'MULTI_RMSD'"
puts $outfile "    cv_ni = 104, cv_nr = 312,"
puts -nonewline $outfile "    cv_i = "
set counter 0
foreach ii [$sel get serial] {
	puts -nonewline $outfile "$ii"
	if {$counter < [llength [$sel get serial]] - 1} {puts -nonewline $outfile ", "}
	incr counter
}
puts $outfile ""
puts -nonewline $outfile "    cv_r = "
set counter 0
foreach x [$sel get x] y [$sel get y] z [$sel get z] {
	puts -nonewline $outfile "    [format "%.3f" $x], [format "%.3f" $y], [format "%.3f" $z]"
	if {$counter < [llength [$sel get serial]] - 1} {puts $outfile ", "}
	incr counter
}
puts $outfile ""
puts $outfile "   npath = 2, path = 5.0, 0, path_mode = 'SPLINE',"
puts $outfile "   nharm = 1, harm = 100"
puts $outfile "/"
close $outfile
quit
EOF

vmd -dispdev text -e tcl > /dev/null
rm tcl
