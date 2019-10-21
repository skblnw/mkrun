#!/bin/bash

TOTAL=20000       # Total number of frames
STEPP=5000

nn=0
rm -f tmpfile
for ii in $(seq 0 $STEPP $TOTAL)
do
    beg=$(( $ii - 10000 ))
    end=$(( $ii - 1     ))
    if [[ $nn -ne 0 ]]
    then
        cat >> tmpfile << EOF
package require pbctools
mol new out.dms
animate delete all
mol addfile run.stk first $beg last $end waitfor all
#pbc unwrap -all
animate write trr $nn.trr
quit
EOF
        vmd -dispdev text -e tmpfile
        rm -f tmpfile
    else
        echo -e "Total number of frame set to $TOTAL. Expect to have $(($TOTAL / $STEPP)) trajectories"
    fi 
    nn=$(($nn + 1))
done

input=""
comma=""
for ii in $(seq 1 $(($TOTAL / $STEPP))); do
    input="$input $ii.trr"
    comma=$comma$'c\n'
done
gmx trjcat -f $input -o out.xtc -settime << EOF
0
$comma
EOF
