#!/bin/bash

LAST=$(qsub job_step5_md-61.pbs)
echo $LAST

for i in {62..70}
do
    LAST=$(qsub -W depend=afterany:$LAST job_step5_md-$i.pbs)
    echo $LAST
done
