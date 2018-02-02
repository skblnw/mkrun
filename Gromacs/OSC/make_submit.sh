#!/bin/bash

LAST=$(qsub -W depend=afterany:2336058 job_step5_md-25.pbs)
echo $LAST

for i in {26..40}
do
    LAST=$(qsub -W depend=afterany:$LAST job_step5_md-$i.pbs)
    echo $LAST
done
