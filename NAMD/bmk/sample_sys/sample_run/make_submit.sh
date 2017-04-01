#!/bin/bash

LAST=$(qsub job-12.pbs)
echo $LAST

for ii in {8,4,3,2,1}
do
    #LAST=$(qsub -W depend=afterany:$LAST job-$ii.pbs)
    LAST=$(qsub job-$ii.pbs)
    echo $LAST
done
