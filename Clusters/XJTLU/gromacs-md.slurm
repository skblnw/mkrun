#!/usr/bin/env bash
#SBATCH --job-name=gmx-md
#SBATCH --account=chunchan
#SBATCH --partition=gpu4090
#SBATCH --qos=4gpus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gpus=1
#SBATCH --mem=32G
#SBATCH --time=6-23:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Inspect the live account/QOS association before copying this example.
set -euo pipefail

tpr=${1:?Usage: sbatch gromacs-md.slurm INPUT.tpr [RUN_PREFIX]}
run_prefix=${2:-${tpr%.tpr}}

[[ -s "$tpr" ]] || { echo "Missing or empty TPR: $tpr" >&2; exit 2; }

module purge
module load gromacs/2023.2-gcc-9.5.0-jzxesel

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Job/account/QOS/partition: ${SLURM_JOB_ID:-unset} ${SLURM_JOB_ACCOUNT:-unset} ${SLURM_JOB_QOS:-unset} ${SLURM_JOB_PARTITION:-unset}"
echo "Visible GPUs: ${CUDA_VISIBLE_DEVICES:-unset}"
module list 2>&1
gmx_mpi --version

mdrun_args=(
  -v
  -s "$tpr"
  -deffnm "$run_prefix"
  -nsteps -1
  -maxh 166.5
  -nb gpu
  -pme gpu
  -bonded gpu
  -update gpu
  -gpu_id 0
)

if [[ -s "${run_prefix}.cpt" ]]; then
  echo "Continuing from ${run_prefix}.cpt"
  mdrun_args+=(-cpi "${run_prefix}.cpt" -append)
fi

srun --ntasks=1 gmx_mpi mdrun "${mdrun_args[@]}"

echo "Finished: $(date --iso-8601=seconds)"
