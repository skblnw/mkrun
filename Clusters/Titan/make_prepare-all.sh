#!/bin/bash

NJOB=4
NN=128
NN_PERJOB=$(($NN / $NJOB))

JOB_DIR=(
    'BAR-PH/tetramer-front/run/'
    'BAR-PH/tetramer-same/run/'
    'BAR-PH/water/run2/'
    'Ups1-Mdm35/system-monomer2/run/'
    )

JOB_INDEX=(103 103 14 310)
    
nheader=28
COMMAND="aprun -n $NN_PERJOB -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0"
for ii in {1..20}; do
    jobname=namd-$ii
    sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
        -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$NN'/g' \
        -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=06:00:00/g' \
        template-job-pbs > job-$ii.pbs
    sed -i $nheader',$d' job-$ii.pbs
    
    echo -e "date\n" >> job-$ii.pbs
    
    for ((jj=0; jj < NJOB; jj++)); do
        echo -e "cd \$PBS_O_WORKDIR/${JOB_DIR[$jj]}" >> job-$ii.pbs
        echo -e "$COMMAND NPT-${JOB_INDEX[$jj]}.namd > NPT-${JOB_INDEX[$jj]}.log &" >> job-$ii.pbs
    done
    
    echo -e "\nwait\ndate" >> job-$ii.pbs
    
    for kk in "${!JOB_INDEX[@]}"; do
        ((JOB_INDEX[$kk]++))
    done
done
