#!/bin/bash

for ii in $(seq 32.0 1.5 54.5)
do
    echo "Extracting for window center $ii"
    awk '{print $1,"",$8}' ../$ii/06_Prod_dist.dat > ../$ii/dist.dat
#    sed -i '$ d' ../$ii/dist.dat
done
