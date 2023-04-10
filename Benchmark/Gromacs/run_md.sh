#!/bin/bash

noupdate=false
gpuid=0
ncore=8
nstep=10000
nn=10
for ii in $ncore; do 
  echo "> Using $ii CPU cores for 1 GPU card and run for $nstep steps"
  if $noupdate; then
    echo -n "> Performing $nn trials of GPU simulations and the benchmark is... "
    total=0
    for jj in $(seq 1 $nn); do
      gmx mdrun -ntomp $ii -gpu_id $gpuid -s md.tpr -deffnm out -nsteps $nstep 2> /dev/null
      res=`grep Performance out.log | awk '{print $2}'`
      total=`echo "${total} + ${res}" | bc`
    done
    echo "scale=1; $total / $nn" | bc
  fi
  echo -n "> Performing $nn trials of GPU (fully offload) simulations and the benchmark is... "
  total=0
  for jj in $(seq 1 $nn); do
    gmx mdrun -ntomp $ii -gpu_id $gpuid -pme gpu -nb gpu -update gpu -s md.tpr -deffnm out -nsteps $nstep 2> /dev/null
    res=`grep Performance out.log | awk '{print $2}'`
    total=`echo "${total} + ${res}" | bc`
  done
  echo "scale=1; $total / $nn" | bc
done
rm \#* out* 
