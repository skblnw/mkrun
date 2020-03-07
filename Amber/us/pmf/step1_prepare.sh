#!/bin/bash

for ii in $(seq 32.0 1.5 60.5)
do
    echo "Extracting for window center $ii"
    awk '{print $1,"",$8}' ../$ii/06_Prod_dist.dat > ../$ii/dist.dat
done
