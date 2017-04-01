#!/bin/bash

#LAST=$(qsub job-16.pbs)
#echo $LAST

for ii in {1,2,4,5,6,7,8,9,10,12}
do
    #LAST=$(qsub -W depend=afterany:$LAST job-$ii.pbs)
    LAST=$(qsub job-$ii.pbs)
    echo $LAST
done
