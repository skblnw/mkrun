#!/bin/bash

for ii in $(seq 32.0 1.5 54.5)
do

echo "Creating histogram for window $ii"
./generate_hist.py -i ../${ii}/dist.dat > hist_${ii}.dat

done
