#!/bin/bash

for ii in 38.0 51.5 53.0 54.5 32.0 33.5 35.0 36.5 56.0 57.5 59.0 60.5
do
    awk '{print $1,"",$8}' ../$ii/06_Prod_dist.dat > ../$ii/dist.dat
done
