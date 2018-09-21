#!/bin/bash

for win in 
do

LAST=$(qsub job-z$win-1.pbs)
echo $LAST

for i in {2..3}
do
    LAST=$(qsub -W depend=afterany:$LAST job-z$win-$i.pbs)
    echo $LAST
done

done
