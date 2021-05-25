#!/bin/bash

mkdir -p result
for ii in 2 4 8 16 32; do
  for jj in $(seq 1 3); do
    gmx mdrun -ntomp $ii -gpu_id 0 -s run_md.tpr -deffnm out -nsteps 10000
    mv out.log result/${ii}-${jj}
  done
  for jj in $(seq 1 3); do
    gmx mdrun -ntomp $ii -pme gpu -nb gpu -update gpu -gpu_id 0 -s run_md.tpr -deffnm out -nsteps 10000
    mv out.log result/${ii}-gpuall-${jj}
  done
done
rm \#* out*
