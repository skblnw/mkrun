#!/bin/bash

NJOB=3
NN=93
NN_PERJOB=$(($NN / $NJOB))

JOB_DIR=(
    'PH-memb/us/2/run/'
    'BAR-PH/water/run1/'
    'BAR-PH/water/run2/'
    )

JOB_INDEX=(1 121 101)
    
COMMAND="aprun -n 12 -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0"
for ii in {1..60}; do
    jobname=namd-$ii
    sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
        -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$NN'/g' \
        -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
        template-job-pbs > job-$ii.pbs
    
    echo -e "date\n" >> job-$ii.pbs
    
    for ((jj=0; jj < NJOB; jj++)); do
        echo -e "cd \$PBS_O_WORKDIR/${JOB_DIR[$jj]}" >> job-$ii.pbs
#        echo -e "$COMMAND NPT-${JOB_INDEX[$jj]}.namd > NPT-${JOB_INDEX[$jj]}.log &" >> job-$ii.pbs
        if [ $jj -eq 0 ]; then
            sed -n '30,/wait/p' ${JOB_DIR[$jj]}/job-us${JOB_INDEX[$jj]}.pbs | head -n -1 >> job-$ii.pbs
        else
            COMMAND="aprun -n 24 -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0"
            echo -e "$COMMAND NPT-${JOB_INDEX[$jj]}.namd > NPT-${JOB_INDEX[$jj]}.log &" >> job-$ii.pbs
        fi
    done
    
    echo -e "\nwait\ndate" >> job-$ii.pbs
    
    for kk in "${!JOB_INDEX[@]}"; do
        ((JOB_INDEX[$kk]++))
    done
done
