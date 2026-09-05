#!/usr/bin/env bash
#SBATCH --job-name=colabfold
#SBATCH --account=chunchan
#SBATCH --partition=gpu4090
#SBATCH --qos=4gpus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

input=${1:?Usage: sbatch colabfold.slurm INPUT.fasta [OUTPUT_DIRECTORY]}
output_dir=${2:-predictions}

[[ -s "$input" ]] || { echo "Missing or empty FASTA: $input" >&2; exit 2; }
mkdir -p -- "$output_dir"

module purge
module load anaconda3
conda_base=$(conda info --base)
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"
conda activate localcolabfold
command -v colabfold_batch >/dev/null || {
  echo "colabfold_batch is not available in the localcolabfold environment." >&2
  exit 2
}

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Job/account/QOS/partition: ${SLURM_JOB_ID:-unset} ${SLURM_JOB_ACCOUNT:-unset} ${SLURM_JOB_QOS:-unset} ${SLURM_JOB_PARTITION:-unset}"
echo "Visible GPUs: ${CUDA_VISIBLE_DEVICES:-unset}"
module list 2>&1
conda info --envs
nvidia-smi

srun --ntasks=1 colabfold_batch "$input" "$output_dir"

echo "Finished: $(date --iso-8601=seconds)"
