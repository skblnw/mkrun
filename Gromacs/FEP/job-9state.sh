#!/bin/bash

# Set some environment variables 
FREE_ENERGY=
echo "Free energy home directory set to $FREE_ENERGY"

MDP=$FREE_ENERGY/MDP
echo ".mdp files are stored in $MDP"

PDB=pdb2gmx
echo "structure files are stored in $FREE_ENERGY/$PDB"

LAMBDA=0

# A new directory will be created for each value of lambda and
# at each step in the workflow for maximum organization.

mkdir Lambda_$LAMBDA
cd Lambda_$LAMBDA

#################################
# ENERGY MINIMIZATION 1: STEEP  #
#################################
echo "Starting minimization for lambda = $LAMBDA..." 

mkdir EM_1 
cd EM_1

# Iterative calls to grompp and mdrun to run the simulations

gmx grompp -f $MDP/EM/em_steep_$LAMBDA.mdp -c $FREE_ENERGY/$PDB/ionized.pdb -p $FREE_ENERGY/$PDB/topol.top -o min$LAMBDA.tpr

gmx mdrun -ntomp 4 -deffnm min$LAMBDA

sleep 10

#####################
# NVT EQUILIBRATION #
#####################
echo "Starting constant volume equilibration..."

cd ../
mkdir NVT
cd NVT

gmx grompp -f $MDP/NVT/nvt_$LAMBDA.mdp -c ../EM_1/min$LAMBDA.gro -p $FREE_ENERGY/$PDB/topol.top -o nvt$LAMBDA.tpr

gmx mdrun -ntomp 4 -deffnm nvt$LAMBDA

echo "Constant volume equilibration complete."

sleep 10

#####################
# NPT EQUILIBRATION #
#####################
echo "Starting constant pressure equilibration..."

cd ../
mkdir NPT
cd NPT

gmx grompp -f $MDP/NPT/npt_$LAMBDA.mdp -c ../NVT/nvt$LAMBDA.gro -p $FREE_ENERGY/$PDB/topol.top -t ../NVT/nvt$LAMBDA.cpt -o npt$LAMBDA.tpr

gmx mdrun -ntomp 4 -deffnm npt$LAMBDA

echo "Constant pressure equilibration complete."

sleep 10

#################
# PRODUCTION MD #
#################
echo "Starting production MD simulation..."

cd ../
mkdir Production_MD
cd Production_MD

gmx grompp -f $MDP/Production_MD/md_$LAMBDA.mdp -c ../NPT/npt$LAMBDA.gro -p $FREE_ENERGY/$PDB/topol.top -t ../NPT/npt$LAMBDA.cpt -o md$LAMBDA.tpr

gmx mdrun -ntomp 4 -deffnm md$LAMBDA

echo "Production MD complete."

# End
echo "Ending. Job completed for lambda = $LAMBDA"
