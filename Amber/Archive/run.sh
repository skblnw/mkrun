#!/bin/bash

export CUDA_VISIBLE_DEVICES=GPUID

pmemd.cuda -O -i 01_min.in -o 01_min.out -p ionized.parm7 -c ionized.rst7 -r 01_min.rst -ref ionized.rst7
pmemd.cuda -O -i 02_heat.in -o 02_heat.out -p ionized.parm7 -c 01_min.rst -r 02_heat.rst -x 02_heat.nc -ref 01_min.rst
pmemd.cuda -O -i 03_eq.in -o 03_eq.out -p ionized.parm7 -c 02_heat.rst -r 03_eq.rst -x 03_eq.nc -ref 02_heat.rst
pmemd.cuda -O -i 04_md.in -o 04_md.out -p ionized.parm7 -c 03_eq.rst -r 04_md.rst -x 04_md.nc -ref 03_eq.rst
#pmemd.cuda -O -i 04_md.in -o PREFIX_md.out -p ionized.parm7 -c 03_eq.rst -r PREFIX_md.rst -x PREFIX_md.nc -ref 03_eq.rst
