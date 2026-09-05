#!/usr/bin/env bash
#SBATCH --job-name=mmpbsa
#SBATCH --account=chunchan
#SBATCH --partition=cpu8358
#SBATCH --qos=52cores
#SBATCH --nodes=1
#SBATCH --ntasks=26
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# This job requests 26 CPU tasks; 52cores is the QOS ceiling, not a requirement
# to occupy all 52 cores. Tune tasks and memory from measured job efficiency.
set -euo pipefail

input=${1:-mmpbsa.in}
complex_tpr=${2:-md.tpr}
index_file=${3:-index.ndx}
trajectory=${4:-md.xtc}
topology=${5:-topol.top}
receptor_group=${REC_GROUP:-1}
ligand_group=${LIG_GROUP:-13}
mmpbsa_env=${MMPBSA_ENV:-gmxmmpbsa}

for path in "$input" "$complex_tpr" "$index_file" "$trajectory" "$topology"; do
  [[ -s "$path" ]] || { echo "Missing or empty input: $path" >&2; exit 2; }
done

module purge
module load anaconda3
conda_base=$(conda info --base)
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"
conda activate "$mmpbsa_env"
command -v gmx_MMPBSA >/dev/null || {
  echo "gmx_MMPBSA is not available in environment: $mmpbsa_env" >&2
  exit 2
}

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Job/account/QOS/partition: ${SLURM_JOB_ID:-unset} ${SLURM_JOB_ACCOUNT:-unset} ${SLURM_JOB_QOS:-unset} ${SLURM_JOB_PARTITION:-unset}"
echo "Tasks: ${SLURM_NTASKS:-unset}"
echo "Receptor/ligand index groups: $receptor_group $ligand_group"
module list 2>&1
conda info --envs

# Index group numbers are system-specific. Confirm them with `gmx make_ndx`
# before submission; the defaults here only make the required choice visible.
srun --ntasks="${SLURM_NTASKS:-26}" gmx_MMPBSA MPI -O \
  -i "$input" \
  -cs "$complex_tpr" \
  -ci "$index_file" \
  -cg "$receptor_group" "$ligand_group" \
  -ct "$trajectory" \
  -cp "$topology"

echo "Finished: $(date --iso-8601=seconds)"
