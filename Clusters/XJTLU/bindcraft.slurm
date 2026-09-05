#!/usr/bin/env bash
#SBATCH --job-name=bindcraft
#SBATCH --account=chunchan
#SBATCH --partition=gpu4090
#SBATCH --qos=4gpus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

settings=${1:?Usage: sbatch bindcraft.slurm SETTINGS.json FILTERS.json ADVANCED.json}
filters=${2:?Usage: sbatch bindcraft.slurm SETTINGS.json FILTERS.json ADVANCED.json}
advanced=${3:?Usage: sbatch bindcraft.slurm SETTINGS.json FILTERS.json ADVANCED.json}
bindcraft_root=${BINDCRAFT_ROOT:-$HOME/opt/BindCraft}

for input in "$settings" "$filters" "$advanced" "$bindcraft_root/bindcraft.py"; do
  [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 2; }
done

module purge
module load anaconda3
conda_base=$(conda info --base)
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"
conda activate BindCraft
command -v python >/dev/null || {
  echo "Python is not available in the BindCraft environment." >&2
  exit 2
}

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Job/account/QOS/partition: ${SLURM_JOB_ID:-unset} ${SLURM_JOB_ACCOUNT:-unset} ${SLURM_JOB_QOS:-unset} ${SLURM_JOB_PARTITION:-unset}"
echo "Visible GPUs: ${CUDA_VISIBLE_DEVICES:-unset}"
echo "BindCraft root: $bindcraft_root"
module list 2>&1
conda info --envs
nvidia-smi

srun --ntasks=1 python -u "$bindcraft_root/bindcraft.py" \
  --settings "$settings" \
  --filters "$filters" \
  --advanced "$advanced"

echo "Finished: $(date --iso-8601=seconds)"
