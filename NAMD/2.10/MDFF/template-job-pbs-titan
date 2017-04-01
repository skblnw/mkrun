#!/bin/bash
#PBS -A bip126
#PBS -N NAME
#PBS -j oe
#PBS -l walltime=02:00:00
#PBS -l nodes=64
##PBS -m e
##PBS -r n
##PBS -W depend=afterany:MM
##
#
#./opt/modules/3.2.6.7/init/bash
source $MODULESHOME/init/bash
#module load fftw
#module load cudatoolkit

#module swap PrgEnv-cray PrgEnv-gnu
module swap PrgEnv-pgi PrgEnv-gnu
module load rca
module load craype-hugepages8M
export HUGETLB_DEFAULT_PAGE_SIZE=8M
export HUGETLB_MORECORE=no

module load namd

export MPICH_PTL_SEND_CREDITS=-1
export MPICH_MAX_SHORT_MSG_SIZE=8000
export MPICH_PTL_UNEX_EVENTS=80000
export MPICH_UNEX_BUFFER_SIZE=100M

cd $PBS_O_WORKDIR
date
aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0 $filename.namd > $filename.log
