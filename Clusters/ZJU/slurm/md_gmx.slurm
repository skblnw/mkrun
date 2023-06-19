#!/bin/bash
#SBATCH -J <job name>
#SBATCH -p multi
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --gres=gpu:1

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"

# Decide the software version
source /public/software/profile.d/apps_gromacs-2020.6.sh

# Assign OMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

MDRUN="gmx mdrun -nb gpu -pme gpu -bonded gpu -update gpu -gpu_id 0"

prefix=t1
$MDRUN -v -s md.tpr -deffnm ${prefix} -nsteps -1 -maxh 72

