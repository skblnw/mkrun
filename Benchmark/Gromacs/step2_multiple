#!/bin/bash

nstep=10000
ncore=8
nn=10
for ii in $ncore; do 
  echo "> Using $ii CPU cores for 1 GPU card and run for $nstep steps"
  if false; then
  echo -n "> Performing $nn trials of GPU simulations and the benchmark is... "
  total=0
  for jj in $(seq 1 $nn); do
    gmx mdrun -ntomp $ii -gpu_id 1 -s run_md.tpr -deffnm out -nsteps $nstep 2> /dev/null
    res=`grep Performance out.log | awk '{print $2}'`
    total=`echo "${total} + ${res}" | bc`
  done
  echo "scale=1; $total / $nn" | bc
  fi
  echo -n "> Performing $nn trials of GPU (fully offload) simulations and the benchmark is... "
  total=0
  for jj in $(seq 1 $nn); do
    gmx mdrun -ntomp $ii -pme gpu -nb gpu -update gpu -gpu_id 1 -s run_md.tpr -deffnm out -nsteps $nstep 2> /dev/null
    res=`grep Performance out.log | awk '{print $2}'`
    total=`echo "${total} + ${res}" | bc`
  done
  echo "scale=1; $total / $nn" | bc
done
rm \#* out* 
