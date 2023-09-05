#!/bin/bash

# Define parameters and commands
GPUID="$1"
GMX="gmx"
MDRUN="gmx mdrun -ntmpi 1 -ntomp 8 -gpu_id $GPUID"
MDRUN_GPU="gmx mdrun -ntmpi 1 -ntomp 8 -gpu_id $GPUID -pme gpu -nb gpu -bonded gpu -update gpu"

# Check if the number of arguments is zero
if [ $# -eq 0 ]; then 
    echo "mkgmx> Usage: $0 [GPU ID]; Suggested: 0"
    exit 1
fi

# Check if the required directory exists
[ ! -d pdb2gmx ] && echo "mkgmx> Directory pdb2gmx not found!" && exit 1 

# Check if the required files exist
files=("pdb2gmx/ionized.gro" "pdb2gmx/topol.top" "pdb2gmx/index.ndx")
for file in "${files[@]}"; do
    [ ! -f "$file" ] && echo "mkgmx> $file not found!" && exit 1
done

# Simulation settings
mkdir -p output
mini=true
double=false
heat=true
eq_npt=true
md1=true
md2=false

INITIAL_PDB=pdb2gmx/ionized.gro
NDX=pdb2gmx/index.ndx
TOP=pdb2gmx/topol.top

# Function to run a simulation step
run_simulation() {
    local previous="$1"
    local prefix="$2"
    local mdp_file="$3"
    local TPR="output/${prefix}.tpr"
    rm -f $TPR
    $GMX grompp -f $mdp_file -o $TPR -c output/${previous}.gro -r $INITIAL_PDB -n $NDX -p $TOP
    $4 -v -s $TPR -deffnm output/${prefix}
}

# Simulation steps based on the set conditions
if $mini; then
    prefix="step1_mini"
    TPR="output/${prefix}.tpr"
    rm -f $TPR
    $GMX grompp -f mdp/step1_mini_double.mdp -o $TPR -c $INITIAL_PDB -r $INITIAL_PDB -n $NDX -p $TOP
    $MDRUN -v -s $TPR -deffnm output/${prefix}
fi
$double && run_simulation "step1_mini" "step1_mini_double" "mdp/step1_mini_double.mdp" "gmx_d mdrun -ntmpi 1 -ntomp 8"
$heat && run_simulation "step1_mini" "step3_annealing" "mdp/step3_annealing.mdp" "$MDRUN"
$eq_npt && run_simulation "step3_annealing" "step4_eq_npt" "mdp/step4_eq_npt.mdp" "$MDRUN"
$md1 && run_simulation "step4_eq_npt" "md" "mdp/step5_md.mdp" "$MDRUN_GPU"

if $md2; then
    previous="md"
    prefix="t1"
    TPR="md.tpr"
    $GMX grompp -f mdp/step5_md.mdp -o $TPR -t output/${previous}.cpt -c $INITIAL_PDB -r $INITIAL_PDB -n $NDX -p $TOP
    # Uncomment the lines below if needed
    # gmx convert-tpr -s ${previous}.tpr -o $prefix.tpr -extend 10000
    # $MDRUN_GPU -v -s $TPR -cpi output/${prefix}.cpt -deffnm output/${prefix} -nsteps -1
    # $MDRUN_GPU -v -s $TPR -deffnm output/${prefix} -nsteps -1
fi
