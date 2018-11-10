#!/bin/bash
#PBS -N NAME
#PBS -j oe
#PBS -l nodes=1:ppn=16

# NAMD binary
export NAMD_PATH=/home/kevin/opt/NAMD_2.10_Linux-x86_64-multicore-CUDA
export LD_LIBRARY_PATH=$NAMD_PATH:$LD_LIBRARY_PATH
export NAMD=$NAMD_PATH/namd2
export CHARMRUN=$NAMD_PATH/charmrun

# CUDA
export CUDA_HOME=/share/apps/cuda
export LD_LIBRARY_PATH=$CUDA_HOME/lib:$LD_LIBRARY_PATH

NP=NUM
NP2=$(expr $NP - 1)



cd $PBS_O_WORKDIR
date
$NAMD +p${NP} +idlepoll bmk-${NP}.namd > bmk-${NP}-v2.log
wait
date