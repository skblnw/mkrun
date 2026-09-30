#!/usr/bin/env bash
#SBATCH --job-name=namd-md
#SBATCH --account=chunchan
#SBATCH --partition=gpu4090
#SBATCH --qos=4gpus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Run one NAMD 3 configuration. The config controls outputName and restarts.
# Usage: sbatch namd-md.slurm INPUT.namd [LOG_PREFIX]
set -euo pipefail

if (( $# < 1 || $# > 2 )); then
  echo "Usage: sbatch namd-md.slurm INPUT.namd [LOG_PREFIX]" >&2
  exit 2
fi

config=$1
log_prefix=${2:-namd-${SLURM_JOB_ID:-manual}}
launcher=/gpfs/work/bio/chunchan/opt/bin/namd3

[[ -s "$config" ]] || { echo "Missing or empty NAMD config: $config" >&2; exit 2; }
[[ -x "$launcher" ]] || { echo "NAMD launcher is unavailable: $launcher" >&2; exit 2; }
[[ ! -e "$log_prefix.log" ]] || { echo "Log already exists: $log_prefix.log" >&2; exit 2; }

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Job/account/QOS/partition: ${SLURM_JOB_ID:-unset} ${SLURM_JOB_ACCOUNT:-unset} ${SLURM_JOB_QOS:-unset} ${SLURM_JOB_PARTITION:-unset}"
echo "Allocated CPUs/visible GPUs: ${SLURM_CPUS_PER_TASK:-unset} ${CUDA_VISIBLE_DEVICES:-unset}"
echo "Launcher: $launcher"
echo "NAMD config/log: $config $log_prefix.log"

# GPU 0 is the first GPU visible inside this one-GPU Slurm allocation.
srun --ntasks=1 "$launcher" +p"${SLURM_CPUS_PER_TASK:-1}" +setcpuaffinity \
  +devices 0 "$config" > "$log_prefix.log" 2>&1

echo "Finished: $(date --iso-8601=seconds)"
