#!/usr/bin/env bash
#SBATCH --job-name=amber-md
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

# Run one Amber stage. The mdin controls new-run versus restart settings.
# Usage: sbatch amber-md.slurm MDIN PARM7 START.rst7 PREFIX [REFERENCE.rst7]
set -euo pipefail

if (( $# < 4 || $# > 5 )); then
  echo "Usage: sbatch amber-md.slurm MDIN PARM7 START.rst7 PREFIX [REFERENCE.rst7]" >&2
  exit 2
fi

mdin=$1
topology=$2
coordinates=$3
prefix=$4
reference=${5:-}
launcher=/gpfs/work/bio/chunchan/opt/bin/pmemd.cuda

for input in "$mdin" "$topology" "$coordinates"; do
  [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 2; }
done
if [[ -n "$reference" ]]; then
  [[ -s "$reference" ]] || { echo "Missing or empty reference: $reference" >&2; exit 2; }
fi
[[ -x "$launcher" ]] || { echo "Amber launcher is unavailable: $launcher" >&2; exit 2; }

# pmemd.cuda -O overwrites files, so refuse a reused prefix.
for output in "$prefix.out" "$prefix.rst7" "$prefix.nc" "$prefix.info"; do
  [[ ! -e "$output" ]] || { echo "Output already exists: $output" >&2; exit 2; }
done

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Job/account/QOS/partition: ${SLURM_JOB_ID:-unset} ${SLURM_JOB_ACCOUNT:-unset} ${SLURM_JOB_QOS:-unset} ${SLURM_JOB_PARTITION:-unset}"
echo "Allocated CPUs/visible GPUs: ${SLURM_CPUS_PER_TASK:-unset} ${CUDA_VISIBLE_DEVICES:-unset}"
echo "Launcher: $launcher"

args=(-O -i "$mdin" -p "$topology" -c "$coordinates"
      -o "$prefix.out" -r "$prefix.rst7" -x "$prefix.nc" -inf "$prefix.info")
if [[ -n "$reference" ]]; then
  args+=(-ref "$reference")
fi
srun --ntasks=1 "$launcher" "${args[@]}"

echo "Finished: $(date --iso-8601=seconds)"
