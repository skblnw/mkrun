#!/bin/bash

LAST=$(qsub job-10.pbs)
echo $LAST

for i in 5 1
do
    LAST=$(qsub -W depend=afterany:$LAST job-$i.pbs)
    echo $LAST
done
