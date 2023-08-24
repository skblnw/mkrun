#!/bin/bash

# Define parameters and commands
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

# Equilibration
equi_prefix="step6.%d_equilibration"
cnt=1
cntmax=6

run_equilibration() {
    local prefix="$1"
    local previous="$2"
    local rest_gro="$3"
    local count="$4"
    
    istep=$(printf ${prefix} ${count})
    pstep=$(printf ${prefix} $((${count}-1)))
    
    if [ ${count} -eq 1 ]; then
        pstep=$mini_prefix
    fi

    # Preprocess the input
    $GMX grompp -f ${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${rest_gro}.gro -p topol.top -n index.ndx

    # Run the equilibration
    $GMX mdrun -v -deffnm ${istep}
}

while [ $cnt -le $cntmax ]; do
    run_equilibration $equi_prefix $mini_prefix $rest_prefix $cnt
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
        $GMX grompp -f ${prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -p topol.top -n index.ndx
    else
        $GMX grompp -f ${prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -t ${pstep}.cpt -p topol.top -n index.ndx
    fi

    # Run the production
    $GMX mdrun -v -deffnm ${istep}
}

while [ $cnt -le $cntmax ]; do
    run_production $prod_step $equi_prefix $prod_prefix $cnt
    cnt=$((cnt+1))
done
