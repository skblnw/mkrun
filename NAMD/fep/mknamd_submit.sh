#!/bin/bash

REPLICA=3

ii=0
for gpu in $(seq 4 7); do
    start=$((ii*REPLICA+1))
    end=$((ii*REPLICA+REPLICA))
    ii=$((ii+1))
    echo "GPU#:$gpu Start:$start End:$end"

    cat > pbs$gpu << EOF
NAMD="namd3 +p1 +devices $gpu"
for ii in \$(seq $start $end); do
    rsync -avh eq/fep.namd eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc trial\$ii
    \$NAMD trial\$ii/fep.namd >& trial\$ii/LOG_fep
done
EOF
    bash pbs$gpu >& LOG$gpu &
done

