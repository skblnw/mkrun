#!/bin/bash
#PBS -N NAME
#PBS -j oe
#PBS -l nodes=12:ppn=16

export IB_PATH=/usr/mpi/gcc/openmpi-1.8.5
export PATH=$IB_PATH/bin:$PATH
export LD_LIBRARY_PATH=$IB_PATH/lib64:$LD_LIBRARY_PATH
MPIRUN=$IB_PATH/bin/mpirun
MPIEXEC=$IB_PATH/bin/mpiexec

#export OPENMPI_PATH=/share/apps/openmpi-1.8.6/bin
#MPIRUN=$OPENMPI_PATH/mpirun
#MPIEXEC=$OPENMPI_PATH/mpiexec

#export NAMD_PATH=/home/kevin/opt/NAMD_2.10_Linux-x86_64-ibverbs-smp-CUDA
#export LD_LIBRARY_PATH=$NAMD_PATH:$LD_LIBRARY_PATH
#export NAMD=$NAMD_PATH/namd2
#export CHARMRUN=$NAMD_PATH/charmrun

export NAMD_PATH=/home/kevin/opt/NAMD_2.10_Linux-x86_64-ibverbs-smp-CUDA
export LD_LIBRARY_PATH=$NAMD_PATH:$LD_LIBRARY_PATH
export NAMD=$NAMD_PATH/namd2
export CHARMRUN=$NAMD_PATH/charmrun

nodefile=/tmp/$PBS_JOBID.nodelist
echo group main ++cpus 16 ++shell ssh > $nodefile
nodes=$( cat $PBS_NODEFILE | awk '!seen[$0]++')
for node in $nodes; do
   echo host $node >> $nodefile
done
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | awk '!seen[$0]++' | wc -l`
NP2=$(expr $NP - $NN)

cd $PBS_O_WORKDIR
date
$CHARMRUN +p${NP2} ++nodelist $nodefile $NAMD cuda.namd ++verbose > cuda_ibverbs_smp_6.log
