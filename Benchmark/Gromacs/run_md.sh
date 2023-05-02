#!/bin/bash

[ $# -eq 0 ] && { echo "> Usage: $0 [# of cores] [GPU ID]"; exit 1; }

noupdate=false
ncore=$1
gpuid=$2
nstep=100000
nn=10

nvidia-smi -L | head -1
lscpu
gmx --version

date
suffix=`date +"%H%M%S"`
for ii in $ncore; do 
  echo "> Using $ii CPU cores for 1 GPU card and run for $nstep steps"
  if $noupdate; then
    echo -n "> Performing $nn trials of GPU simulations and the benchmark is... "
    total=0
    for jj in $(seq 1 $nn); do
      gmx mdrun -ntomp $ii -gpu_id $gpuid -s md.tpr -deffnm out-${suffix} -nsteps $nstep 2> /dev/null
      res=`grep Performance out-${suffix}.log | awk '{print $2}'`
      echo $res >> nsday-${suffix}
      total=`echo "${total} + ${res}" | bc`
    done
    echo "scale=1; $total / $nn" | bc
  fi
  echo -n "> Performing $nn trials of GPU (fully offload) simulations and the benchmark is... "
  total=0
  for jj in $(seq 1 $nn); do
    gmx mdrun -ntomp $ii -gpu_id $gpuid -pme gpu -nb gpu -update gpu -s md.tpr -deffnm out-${suffix} -nsteps $nstep 2> /dev/null
    res=`grep Performance out-${suffix}.log | awk '{print $2}'`
    echo $res >> nsday-${suffix}
    total=`echo "${total} + ${res}" | bc`
  done
  echo "scale=1; $total / $nn" | bc
done

date
rm \#* out-${suffix}* 
