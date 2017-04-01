#!/bin/bash
#PBS -N mart
#PBS -j oe
#PBS -l nodes=1:ppn=16

export OPENMPI_PATH=/share/apps/openmpi-1.8.6/bin
MPIRUN=$OPENMPI_PATH/mpirun
MPIEXEC=$OPENMPI_PATH/mpiexec

#export GMX_PATH=/share/apps/gromacs-5.1.1-MPI-double/
export GMX_PATH=/home/kevin/opt/gromacs-5.1.1-MPI-CUDA-single/
source $GMX_PATH/bin/GMXRC.bash
GMXBIN=$GMX_PATH/bin/gmx_mpi

NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | awk '!seen[$0]++' | wc -l`
NP2=$(expr $NP - $NN)

# Number of Proccessors Per Node
NPPN=16

cd $PBS_O_WORKDIR
date

mini=false
eq=false
md1=false
md2=false

cd output

if $mini; then
    for ii in {1..1}; do
        jj=$(expr $ii - 1)
        if [ "$ii" -eq 1 ]; then
            $GMXBIN grompp -f ../mini-steep.mdp -c ../../system.gro -p ../../system.top -n ../index.ndx -o mini-$ii.tpr
            wait
        else 
            $GMXBIN grompp -f ../mini-steep.mdp -c mini-$jj.gro -p ../../system.top -n ../index.ndx -o mini-$ii.tpr
            wait
        fi
        $GMXBIN mdrun -v -ntomp $NPPN -pin on -deffnm mini-$ii
        wait
    done
fi

if $eq; then
    mini=1
    $GMXBIN grompp -f ../eq.mdp -c mini-$mini.gro -p ../../system.top -n ../index.ndx -o eq.tpr
    wait
    $GMXBIN mdrun -v -ntomp $NPPN -pin on -deffnm eq
    wait
fi

if $md1; then
    $GMXBIN grompp -f ../md.mdp -c eq.gro -p ../../system.top -n ../index.ndx -o md-1.tpr
    wait
    $GMXBIN mdrun -v -ntomp $NPPN -pin on -deffnm md-1
    wait
fi

if $md2; then
    for ii in {2..2}
    do
        jj=$(expr $ii - 1)
        $GMXBIN convert-tpr -s md-$jj.tpr -o md-$ii.tpr -extend 500000
        wait
        $GMXBIN mdrun -v -ntomp $NPPN -pin on -s md-$ii.tpr -cpi md-$jj.cpt -deffnm md-$ii
        wait
    done
fi

date
