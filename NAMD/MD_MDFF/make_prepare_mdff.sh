#!/bin/bash
################################################
# Aug 2015
# CBMSG
# Kevin
# This script was created for preparing NAMD MD jobs ***IN A LAZY WAY***
#
# Updates:
# May 2017
#  Revised NAMD CUDA-smp-ib commands (highly recommended in the release note)
#  Deleted NVT equilibirum before releasing constraints
#  Added releasing constraints to lipid heads
#  Deleted redundant minimization steps
# 17 Aug 2016
#  Created mknamd_prepare
#  Added function make_pbs
#  Added clusters as options
# 13 Sep 2015
#  Added option bmk
# 4 Aug 2015
#  Added variable previous
# 3 Aug 2015
#  Added variable COMMAND
#  Added variable jobname
#  Added variable nnode and nheader
################################################
# NAMD Commands (GPU):
#   COMBO (remote-shell=ssh): \$CHARMRUN +p\$NP_2 ++ppn 15 ++scalable-start ++nodelist \$nodefile \$NAMD +setcpuaffinity +pemap 1-15 +commap 0 ++verbose
#   COMBO (remote-shell=ib): \$CHARMRUN +p\$NP_2 ++ppn 15 ++scalable-start ++verbose ++mpiexec ++remote-shell ./mpiexec_wrapper \$NAMD +idlepoll +setcpuaffinity +pemap 1-15 +commap 0
#   Titan: aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0
#   College: mpirun -n \$NP \$NAMD
################################################
# NAMD Commands (CPU):
#   COMBO: \$CHARMRUN +p\$NP ++ppn 16 ++scalable-start ++nodelist \$nodefile \$NAMD ++verbose
################################################
# Gromacs Commands (GPU):
#   COMBO: $MPIRUN -mca btl self,openib -np $NN -npernode 1 $GMXBIN mdrun -v -ntomp 7 -pin on -s W50k.tpr -deffnm output/gpu-1node-NNomp
################################################

mkdir -p output
mkdir -p log
mkdir -p restraints

nnode=8
ncpu=10
ncpu_1=$(($ncpu-1))
ncpu_2=$(($ncpu-2))

# Cluster Choosing
# true if you choose that cluster
combo=true
# Combo has a cpu version of pbs
combo_cpu=false
titan=false
# Titan has a multiple version of pbs
titan_multiple=false
college=false
if $combo; then
    if [[ $nnode -eq 1 ]]; then
        COMMAND="\$NAMD +p$ncpu_1"
    elif $combo_cpu; then
        COMMAND="\$CHARMRUN +p$ncpu_1 ++scalable-start ++verbose ++mpiexec ++remote-shell ./mpiexec_wrapper \$NAMD"
    else
        COMMAND="\$CHARMRUN +p\${NP_2} ++ppn $ncpu_2 ++scalable-start ++verbose ++mpiexec ++remote-shell ./mpiexec_wrapper \$NAMD +idlepoll +setcpuaffinity +pemap 1-$ncpu_2 +commap 0"
    fi
elif $titan; then
    COMMAND="aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0"
elif $college; then
    COMMAND="mpirun -n \$NP \$NAMD"
else
    echo "Specify on which machine you are gonna run!"
    exit 0
fi


membrane_exist=true
# MD Steps Choosing
# true to turn on
mdff=true

# Check existence of a file
check_exist () {
    if [ ! -s "$1" ]; then
        echo "$1 does not exist! Be careful!"
        exit 1
    fi
}

# Create a pbs acoording to chosen cluster
# Usage: make_pbs $jobname ($nnode_total)
make_pbs () {
    check_exist template-job-pbs

    if [ -z "$1" ]; then
        echo "-Parameter \$jobname is empty!"
        exit 1
    fi

    sed -e 's/^#PBS -N .*/#PBS -N '$1'/g' template-job-pbs > job-$1.pbs

    if $titan; then
        if [ -z "$2" ]; then
            sed -i 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode'/g' job-$1.pbs
        else
            sed -i 's/^#PBS -l nodes=.*/#PBS -l nodes='$2'/g' job-$1.pbs
        fi
    else
        sed -i 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' job-$1.pbs
    fi
}


if $mdff; then
    for ii in {1..1}
    do
        jobname=mdff$ii
        make_pbs $jobname
        outputname=symmetry

        jj=$(expr $ii - 1)
        echo -e "$COMMAND ${outputname}-step$ii.namd > log/${outputname}-step$ii.log\nwait\ndate" >> job-$jobname.pbs
    done
fi
