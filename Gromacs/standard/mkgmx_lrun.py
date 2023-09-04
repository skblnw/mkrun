#!/usr/bin/env python3

import argparse
import os
import subprocess

def check_file_existence(filename):
    """Checks if the specified file exists."""
    if not os.path.isfile(filename):
        print(f"mkvmd> {filename} not found!")
        exit(1)

def execute_subprocess(cmd):
    """Executes the provided command using subprocess."""
    subprocess.run(cmd, shell=True)

def start_simulation(prefix, tpr_file, mdrun_gpu, nsteps):
    """Start a new molecular dynamics simulation."""
    check_file_existence(f"{prefix}.cpt")
    print("mkvmd> Starting a new simulation")
    cmd = f"{mdrun_gpu} -v -s {tpr_file} -deffnm {prefix} -nsteps {nsteps}"
    execute_subprocess(cmd)

def continue_simulation(prefix, previous, tpr_file, mdrun_gpu, nsteps):
    """Continue a previous molecular dynamics simulation."""
    check_file_existence(f"{previous}.cpt")
    # Backing up the previous checkpoint file
    os.system(f"cp {previous}.cpt {previous}.cpt.BAK")
    print(f"mkvmd> Continuing simulation from {previous}.cpt")
    cmd = f"{mdrun_gpu} -v -s {tpr_file} -cpi {previous}.cpt -deffnm {prefix} -nsteps {nsteps}"
    execute_subprocess(cmd)

def main():
    parser = argparse.ArgumentParser(description='Start or continue a molecular dynamics simulation.')
    parser.add_argument('gpuid', type=int, help='GPU ID to use')
    parser.add_argument('prefix', type=str, help='Prefix for the current run')
    parser.add_argument('--previous', type=str, default="", help='Prefix for the previous run (only for continuation)')
    parser.add_argument('--tpr', type=str, default="md.tpr", help='Name of the TPR file (default is "md.tpr")')
    parser.add_argument('--nsteps', type=int, default=-1, help='Number of steps for the simulation (default is -1)')
    
    args = parser.parse_args()

    MDRUN_GPU = f"gmx mdrun -ntmpi 1 -ntomp 16 -gpu_id {args.gpuid} -pme gpu -nb gpu -bonded gpu -update gpu"

    # Check if the TPR file exists
    check_file_existence(args.tpr)

    # Start a new simulation or continue from previous checkpoint
    if not args.previous:
        start_simulation(args.prefix, args.tpr, MDRUN_GPU, args.nsteps)
    else:
        continue_simulation(args.prefix, args.previous, args.tpr, MDRUN_GPU, args.nsteps)

if __name__ == "__main__":
    main()
