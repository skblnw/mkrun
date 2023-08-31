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
files=("pdb2gmx/ionized.pdb" "pdb2gmx/topol.top" "pdb2gmx/index.ndx")
for file in "${files[@]}"; do
    [ ! -f "$file" ] && echo "mkgmx> $file not found!" && exit 1
done

# Simulation settings
mini=false
eq=true
md=false

TOP="pdb2gmx/topol.top"
NDX="pdb2gmx/index.ndx"
INITIAL_PDB="pdb2gmx/ionized.pdb"

# Minimization
mini_prefix="step6.0_minimization"

run_minimization() {
    local prefix="$1"
    local init="$2"
    
    # Preprocess the input
    $GMX grompp -f mdp/${prefix}.mdp -o ${prefix}.tpr -c ${INITIAL_PDB} -r ${INITIAL_PDB} -p ${TOP} -n ${NDX}

    # Run the minimization
    $MDRUN -v -deffnm ${prefix}
}

$mini && run_minimization $mini_prefix

# Equilibration
equi_prefix="step6.%d_equilibration"
cnt=1
cntmax=6

run_equilibration() {
    local prefix="$1"
    local count="$2"
    
    istep=$(printf ${prefix} ${count})
    pstep=$(printf ${prefix} $((${count}-1)))
    
    if [ ${count} -eq 1 ]; then
        pstep=$mini_prefix
    fi

    # Preprocess the input
    $GMX grompp -f mdp/${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${INITIAL_PDB} -p ${TOP} -n ${NDX}
    # Run the equilibration
    $MDRUN -v -deffnm ${istep}
}

while [ $cnt -le $cntmax ]; do
    $eq && run_equilibration $equi_prefix $cnt
    cnt=$((cnt+1))
done

# Production
prod_step="step7"
prod_prefix="step7_production"
cnt=1
cntmax=1

run_production() {
    local step="$1"
    local previous="$2"
    local prefix="$3"
    local count="$4"
    
    istep="${step}_${count}"
    pstep="${step}_$((${count}-1))"
    
    if [ ${count} -eq 1 ]; then
        pstep=$(printf ${equi_prefix} 6)
        $GMX grompp -f ${prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -p ${TOP} -n ${NDX}
    else
        $GMX grompp -f ${prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -t ${pstep}.cpt -p ${TOP} -n ${NDX}
    fi

    # Run the production
    $MDRUN_GPU -v -deffnm ${istep}
}

while [ $cnt -le $cntmax ]; do
    $md && run_production $prod_step $equi_prefix $prod_prefix $cnt
    cnt=$((cnt+1))
done
