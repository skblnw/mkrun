#!/bin/bash
#PBS -N PREFIX-NN
#PBS -A PAS1326
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=28:gpus=1
##PBS -l mem=120GB
#PBS -j oe
##PBS -q parallel

module load gromacs
export OMP_NUM_THREADS=24
export MV2_ENABLE_AFFINITY=0
PPN=1
NTOMP=$OMP_NUM_THREADS

cd $PBS_O_WORKDIR

prefix=PREFIX
iname=INPUTNAME
oname=OUTPUTNAME
CONTINUE=false
if $CONTINUE; then
#    gmx_mpi grompp -f $prefix.mdp -o $oname.tpr -c pdb2gmx/ionized.gro -t output/$iname.cpt -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    gmx convert-tpr -s $iname.tpr -o $oname.tpr -extend 5000
    mpiexec -ppn $PPN gmx_gpu_mpi mdrun -ntomp $NTOMP -pin on -gpu_id 0 -v -s $oname.tpr -cpi output/$iname.cpt -deffnm output/$oname
else
    gmx_mpi grompp -f $prefix.mdp -o $oname.tpr -c output/$iname.gro -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    mpiexec -ppn $PPN gmx_gpu_mpi mdrun -ntomp $NTOMP -pin on -gpu_id 0 -v -s $oname.tpr -deffnm output/$oname
fi
