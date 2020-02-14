#!/bin/bash

export prmtop=ionized.parm7
pmemd.cuda -O -i 01_min.in -o 01_min.out -p $prmtop -c ionized.rst7 -r 01_min.rst -ref ionized.rst7
pmemd.cuda -O -i 02_heat.in -o 02_heat.out -p $prmtop -c 01_min.rst -r 02_heat.rst -x 02_heat.nc -ref 01_min.rst
pmemd.cuda -O -i 03_eq.in -o 03_eq.out -p $prmtop -c 02_heat.rst -r 03_eq.rst -x 03_eq.nc -ref 02_heat.rst
pmemd.cuda -O -i 06_Prod.in -o 06_Prod.out -p $prmtop -c 03_eq.rst -r 06_Prod.rst -x 06_Prod.nc
