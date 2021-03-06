#!/bin/bash

TOTAL=100000       # Total number of frames

nn=0
rm -f tmpfile
for ii in $(seq 0 10000 $TOTAL)
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
        echo -e "Total number of frame set to $TOTAL. Expect to have $(($TOTAL / 10000)) trajectories"
    fi 
    nn=$(($nn + 1))
done
