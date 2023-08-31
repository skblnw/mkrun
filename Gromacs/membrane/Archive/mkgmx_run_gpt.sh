#!/bin/bash

# Define parameters and commands
<<<<<<< HEAD
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
=======
GMX="gmx"

# Minimization
mini_prefix="step6.0_minimization"
init="step5_input"
rest_prefix="step5_input"

run_minimization() {
    local prefix="$1"
    local init_gro="$2"
    local rest_gro="$3"
    
    # Preprocess the input
    $GMX grompp -f ${prefix}.mdp -o ${prefix}.tpr -c ${init_gro}.gro -r ${rest_gro}.gro -p topol.top -n index.ndx

    # Run the minimization
    gmx_d mdrun -v -deffnm ${prefix}
}

run_minimization $mini_prefix $init $rest_prefix
>>>>>>> 6969ea4c94249aa83711e1a83638b393307ad7e0

# Equilibration
equi_prefix="step6.%d_equilibration"
cnt=1
cntmax=6

run_equilibration() {
    local prefix="$1"
<<<<<<< HEAD
    local count="$2"
=======
    local previous="$2"
    local rest_gro="$3"
    local count="$4"
>>>>>>> 6969ea4c94249aa83711e1a83638b393307ad7e0
    
    istep=$(printf ${prefix} ${count})
    pstep=$(printf ${prefix} $((${count}-1)))
    
    if [ ${count} -eq 1 ]; then
        pstep=$mini_prefix
    fi

    # Preprocess the input
<<<<<<< HEAD
    $GMX grompp -f mdp/${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${INITIAL_PDB} -p ${TOP} -n ${NDX}
    # Run the equilibration
    $MDRUN -v -deffnm ${istep}
}

while [ $cnt -le $cntmax ]; do
    $eq && run_equilibration $equi_prefix $cnt
=======
    $GMX grompp -f ${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${rest_gro}.gro -p topol.top -n index.ndx

    # Run the equilibration
    $GMX mdrun -v -deffnm ${istep}
}

while [ $cnt -le $cntmax ]; do
    run_equilibration $equi_prefix $mini_prefix $rest_prefix $cnt
>>>>>>> 6969ea4c94249aa83711e1a83638b393307ad7e0
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
<<<<<<< HEAD
        $GMX grompp -f ${prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -p ${TOP} -n ${NDX}
    else
        $GMX grompp -f ${prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -t ${pstep}.cpt -p ${TOP} -n ${NDX}
    fi

    # Run the production
    $MDRUN_GPU -v -deffnm ${istep}
}

while [ $cnt -le $cntmax ]; do
    $md && run_production $prod_step $equi_prefix $prod_prefix $cnt
=======
        $GMX grompp -f ${prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -p topol.top -n index.ndx
    else
        $GMX grompp -f ${prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -t ${pstep}.cpt -p topol.top -n index.ndx
    fi

    # Run the production
    $GMX mdrun -v -deffnm ${istep}
}

while [ $cnt -le $cntmax ]; do
    run_production $prod_step $equi_prefix $prod_prefix $cnt
>>>>>>> 6969ea4c94249aa83711e1a83638b393307ad7e0
    cnt=$((cnt+1))
done
