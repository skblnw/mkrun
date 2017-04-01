#!/bin/bash

LAST=$(qsub job-npt2.pbs)
echo $LAST

for i in {3..5}
do
    LAST=$(qsub -W depend=afterany:$LAST job-npt$i.pbs)
    echo $LAST
done
