#!/usr/bin/env bash
#SBATCH --job-name=charmmgui
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

# MODE=equilibration runs the CHARMM-GUI README stages once.
# MODE=production uses STEP (step5, step7, ...) and RUN_PREFIX.
set -euo pipefail

mode=${MODE:-equilibration}
step=${STEP:-step7}
run_prefix=${RUN_PREFIX:-t1}

module purge
module load gromacs/2023.2-gcc-9.5.0-jzxesel

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Job/account/QOS/partition: ${SLURM_JOB_ID:-unset} ${SLURM_JOB_ACCOUNT:-unset} ${SLURM_JOB_QOS:-unset} ${SLURM_JOB_PARTITION:-unset}"
echo "Visible GPUs: ${CUDA_VISIBLE_DEVICES:-unset}"
echo "Mode/step/prefix: $mode $step $run_prefix"
module list 2>&1
gmx_mpi --version

case "$mode" in
  equilibration)
    [[ -s README ]] || {
      echo "CHARMM-GUI's generated README is missing or empty." >&2
      exit 2
    }
    command -v csh >/dev/null || {
      echo "csh is required for the generated CHARMM-GUI commands." >&2
      exit 2
    }

    wrapper="charmmgui.${SLURM_JOB_ID:-manual}.csh"
    trap 'rm -f -- "$wrapper"' EXIT

    # Keep the generated instructions intact; select GROMACS MPI, request GPU
    # offload for mdrun, and run each equilibration stage once. Inspect the
    # wrapper before adapting this to a differently generated README.
    sed -e 's/gmx mdrun -v/gmx_mpi mdrun -gpu_id 0 -nb gpu -pme gpu -bonded gpu -v/g' \
        -e 's/gmx_d mdrun /gmx_mpi mdrun /g' \
        -e 's/gmx grompp /gmx_mpi grompp /g' \
        -e 's/^set cntmax = 10/set cntmax = 1/' \
        README > "$wrapper"
    csh "$wrapper"
    ;;

  production)
    mdp="${step}_production.mdp"
    start_gro="${step}_1.gro"
    start_cpt="${step}_1.cpt"
    for input in "$mdp" "$start_gro" "$start_cpt" index.ndx topol.top; do
      [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 2; }
    done

    if [[ ! -s md.tpr ]]; then
      gmx_mpi grompp \
        -f "$mdp" \
        -c "$start_gro" \
        -t "$start_cpt" \
        -n index.ndx \
        -p topol.top \
        -o md.tpr
    fi

    mdrun_args=(
      -v
      -s md.tpr
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
    ;;

  *)
    echo "MODE must be equilibration or production, not: $mode" >&2
    exit 2
    ;;
esac

echo "Finished: $(date --iso-8601=seconds)"
