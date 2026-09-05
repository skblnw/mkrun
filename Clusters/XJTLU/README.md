# XJTLU Slurm examples

These files are compact examples for manual inspection. They document working
patterns observed on the XJTLU cluster; they are not universal submission
wrappers. The partition, QOS, account, module, environment and input choices
must be checked before each real submission.

Cluster snapshot: September 2026.

## Examples

- `gromacs-md.slurm`: a checkpoint-aware, single-GPU GROMACS production run.
- `gromacs-charmmgui.slurm`: CHARMM-GUI equilibration or production without
  maintaining separate scripts for every numbered step or trajectory format.
- `colabfold.slurm`: a single-GPU LocalColabFold prediction.
- `bindcraft.slurm`: a single-GPU BindCraft design run.
- `mmpbsa-cpu.slurm`: a CPU/MPI `gmx_MMPBSA` analysis.

The examples use accounting account `chunchan`. Change `--account` to the
association that should be charged, including any account-name variation such
as a project-specific suffix. An account, a QOS and a resource request are
different things:

- `--account` chooses where usage is accounted.
- `--qos=4gpus` selects an entitlement/concurrency policy. It does **not** make
  one job consume four GPUs.
- `--gpus=1` is the GPU request made by that job.

## Current cluster shape

At the time of inspection:

| Workload | Partitions | Representative production QOS |
| --- | --- | --- |
| RTX GPU | `gpu3090`, `gpu4090` | `1gpu`, `2gpu`, `4gpus`, `6gpus`, `8gpus`, `12gpus`, `16gpus` |
| A800 GPU | `gpua800` | `1a800`, `2a800`, `3a800`, `4a800`, `8a800`, `16a800` |
| CPU | `cpu6348`, `cpu8358` | for example `52cores` |
| Large-memory CPU | `cpuhuge1t` | for example `46cores` |

Do not copy an A800 choice unless the selected account is actually associated
with the matching A800 QOS. Debug QOS entries are intentionally omitted here.
Cluster policy changes, so verify the live state rather than treating this
table as permanent:

```bash
sacctmgr show assoc where user="$USER"
scontrol show partition -o
module spider gromacs
```

## Principles used in these examples

1. State the account, partition, QOS, wall time, memory, logs and requested
   CPUs/GPUs explicitly. A QOS ceiling is not a reason to request its maximum.
2. Give jobs meaningful names and keep `%x` (job name) and `%j` (job ID) in
   log filenames so concurrent runs do not overwrite each other.
3. Fail early when inputs, commands or environments are missing. Print the
   host, allocation and software version into the log for provenance.
4. Resume long simulations from checkpoints. For multi-stage work, prefer
   `sbatch --dependency=afterok:<jobid>` so a failed stage does not start the
   next one. Use arrays for genuinely independent replicas or parameter sets.
5. Inspect completed jobs instead of tuning from requested resources alone:

   ```bash
   sacct -j JOBID --format=JobID,State,Elapsed,AllocTRES,MaxRSS,ExitCode
   seff JOBID
   ```

6. Do not assume a scratch directory or copy-back policy until the cluster's
   supported scratch location and retention rules have been confirmed.

## Example submissions

```bash
sbatch gromacs-md.slurm md.tpr replica-01

sbatch --export=ALL,MODE=equilibration \
  gromacs-charmmgui.slurm

sbatch --export=ALL,MODE=production,STEP=step7,RUN_PREFIX=replica-01 \
  gromacs-charmmgui.slurm

sbatch colabfold.slurm targets.fasta predictions

sbatch bindcraft.slurm settings.json filters.json advanced.json

sbatch --export=ALL,REC_GROUP=1,LIG_GROUP=13 \
  mmpbsa-cpu.slurm mmpbsa.in md.tpr index.ndx md.xtc topol.top
```

The CHARMM-GUI example replaces the previous step-specific and output-specific
duplicates. Trajectory settings belong in the selected MDP file rather than in
another scheduler script. The old standalone `trjconv` helper was also removed
because this directory is now limited to Slurm examples and their explanation.
