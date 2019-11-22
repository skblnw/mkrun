#!/bin/bash

TOTAL=10000       # Total number of frames
STEPP=1000        # Number of frames for each trajectory fragment
TSTEP=120         # Timestep in ps

nn=0
rm -f tmpfile
for ii in $(seq 0 $STEPP $TOTAL)
do
    beg=$(( $ii - $STEPP ))
    end=$(( $ii - 1     ))
    if [[ $nn -ne 0 ]]
    then
        echo $nn $beg $end
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
for ii in $(seq 1 $(($TOTAL / $STEPP))); do
    input="$input $ii.trr"
done
gmx trjcat -f $input -o out.xtc -cat
gmx trjconv -f out.xtc -o pf${TSTEP}ps.xtc -timestep $TSTEP