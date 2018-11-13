#!/bin/bash
# gmx mdrun -v -deffnm md -nsteps 10000

for ii in 2 8 16 32
do

  sed -e 's/NN/'$ii'/g' template-pbs > bmk.pbs
  qsub bmk.pbs

done