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

# Number of Proccessors Per Node
NP=16

COMMAND="mpirun -n $NP -hostfile $PBS_NODEFILE $GMXBIN"

cd $PBS_O_WORKDIR
date

flexible=false
mini=false
nvt=false
npt=false
md1=false
md2=false

cd output

if $flexible; then
    gmx_mpi grompp -c ../../system.gro -f ../bmw_mini.mdp -p ../../flexible -n ../index -o flexible

    $COMMAND mdrun -table ../table -deffnm flexible
fi

if $mini; then
    for ii in {1..1}; do
        jj=$(expr $ii - 1)
        if [ "$ii" -eq 1 ]; then
            $GMXBIN grompp -f ../bmw_mini.mdp -c flexible.gro -p ../../system -n ../index.ndx -o mini-$ii.tpr
            wait
        else 
            $GMXBIN grompp -f ../bmw_mini.mdp -c mini-$jj.gro -p ../../system -n ../index.ndx -o mini-$ii.tpr
            wait
        fi
        $COMMAND mdrun -v -table ../table -deffnm mini-$ii
        wait
    done
fi

if $nvt; then
    mini=1
    $GMXBIN grompp -f ../bmw_nvt.mdp -c mini-1.gro -p ../../system -n ../index.ndx -o nvt.tpr
    wait
    $COMMAND mdrun -v -table ../table -deffnm nvt
    wait
fi

if $npt; then
    mini=1
    $GMXBIN grompp -f ../bmw_npt.mdp -c mini-1.gro -p ../../system -n ../index.ndx -o npt.tpr
    wait
    $COMMAND mdrun -v -table ../table -deffnm npt
    wait
fi

if $md1; then
    $GMXBIN grompp -f ../bmw_md.mdp -c npt.gro -p ../../system -n ../index.ndx -o md-1.tpr
    wait
    $COMMAND mdrun -v -cpi npt.cpt -nsteps 5000 -table ../table -deffnm md-1
    wait
fi

if $md2; then
    for ii in {STARTJ..ENDJ}
    do
        jj=$(expr $ii - 1)
        $GMXBIN convert-tpr -s md-$jj.tpr -o md-$ii.tpr -extend 5000000
        wait
        $COMMAND mdrun -v -table ../table -s md-$ii.tpr -cpi md-$jj.cpt -deffnm md-$ii
        wait
    done
fi

date
