#!/bin/bash

perform_benchmark() {
    local n_trials=$1
    local n_cores=$2
    local devices=$3
    local log_file=$4
    local grep_pattern=$5
    local awk_field=$6
    local namd_file=$7

    local total=0
    for trial in $(seq 1 $n_trials); do
        namd3 +p$n_cores +devices $devices $namd_file > $log_file
        local result=$(grep $grep_pattern $log_file | awk -v field=$awk_field '{sum+=$field} END {print sum/NR}')
        total=$(echo "${total} + ${result}" | bc)
    done

    echo "scale=1; $total / $n_trials" | bc
}

cpu_cores=1
num_trials=10

for cpu_core in $cpu_cores; do 
    
    log_file="LOG"
    
    if [ $cpu_core -ne 1 ]; then
        echo "> Using $cpu_core CPU cores for 1 GPU card"
        echo -n "> Performing $n_trials trials of GPU simulations and the benchmark is... "

        perform_benchmark $num_trials $cpu_cores 0 $log_file "Benchmark" 8 "fep.soft.namd"
    fi
    
    echo -n "> Performing $num_trials trials of GPU (Single-Node, fully offload) simulations and the benchmark is... "
    
    perform_benchmark $num_trials 1 0 $log_file "TIMING:" 9 "fep.soft.gpu.namd"

    rm $log_file
done

rm FFTW_* alchemy*
