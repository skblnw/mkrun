#!/bin/bash
alias vmd='rlwrap /Applications/VMD\ 1.9.4a38.app/Contents/MacOs/startup.command'

for ii in $(seq 1 32)
do
    cd $ii
    cat > cvin.tcl << 'EOF'
mol new ionized.parm7.NN type parm7
mol addfile ionized.rst7.NN type rst7
set outf [open "cv.in" "w"]
puts $outf "&colvar cv_type='MULTI_RMSD',"

set sel [atomselect top "protein and resid 892 to 1078 1254 to 1300 and name CA"]
puts $outf "cv_ni=[$sel num],"
puts $outf "cv_nr=[expr [$sel num] * 3],"
puts -nonewline $outf "cv_i="
foreach ii [$sel get index] {
    set aa [atomselect top "index $ii"]
    puts -nonewline $outf "[$aa get serial],"
}
puts -nonewline $outf "\ncv_r="
foreach ii [$sel get index] {
    set aa [atomselect top "index $ii"]
    puts $outf "[$aa get x], [$aa get y], [$aa get z],"
}
puts $outf "\nanchor_position=-3,0,0,3,\nanchor_strength=5,5, /"
close $outf
quit
EOF
    sed -i.BAK "s/NN/$ii/g" cvin.tcl
    vmd -dispdev text -e cvin.tcl
    rm -f cvin.tcl cvin.tcl.BAK
    cd ..
done
