#!/bin/bash
##SBATCH -J <job name>
#SBATCH -p single
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"

# Decide the software version
source /public/software/profile.d/apps_namd-3.0alpha9.sh

NAMD="/public/software/apps/NAMD_3.0alpha9/namd3"

$NAMD +p1 +devices 0 eq/fep.eq.namd >& eq/LOG_eq

nn=1
rsync -avh eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc eq/fep.namd trial${nn}
# 24 windows
sed -i 's/^set all.*/set all { 0 0.00001 0.0001 0.001 0.01 0.02 0.04 0.06 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.94 0.96 0.98 0.99 0.999 0.9999 0.99999 1.00 }/g' trial${nn}/fep.namd
sed -i 's/^set numSteps.*/set numSteps 600000/g'  trial${nn}/fep.namd
$NAMD +p1 +devices 0 trial${nn}/fep.namd >& trial${nn}/LOG_fep

nn=2
rsync -avh eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc eq/fep.namd trial${nn}
# 38 windows
sed -i 's/^set all.*/set all { 0 0.00001 0.0001 0.001 0.004 0.01 0.02 0.04 0.06 0.1 0.14 0.18 0.22 0.26 0.30 0.34 0.38 0.42 0.46 0.5 0.54 0.58 0.62 0.66 0.7 0.74 0.78 0.82 0.86 0.9 0.94 0.96 0.98 0.99 0.996 0.999 0.9999 0.99999 1.00 }/g' trial${nn}/fep.namd
sed -i 's/^set numSteps.*/set numSteps 450000/g'  trial${nn}/fep.namd
$NAMD +p1 +devices 0 trial${nn}/fep.namd >& trial${nn}/LOG_fep

nn=3
rsync -avh eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc eq/fep.namd trial${nn}
# 60 windows
sed -i 's/^set all.*/set all { 0 0.000001 0.00001 0.0001 0.001 0.002 0.004 0.006 0.008 0.01 0.014 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.12 0.14 0.18 0.22 0.26 0.30 0.34 0.38 0.42 0.46 0.5 0.54 0.58 0.62 0.66 0.7 0.74 0.78 0.82 0.86 0.88 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.986 0.99 0.992 0.994 0.996 0.998 0.999 0.9999 0.99999 0.999999 1.00 }/g' trial${nn}/fep.namd
sed -i 's/^set numSteps.*/set numSteps 300000/g'  trial${nn}/fep.namd
$NAMD +p1 +devices 0 trial${nn}/fep.namd >& trial${nn}/LOG_fep
