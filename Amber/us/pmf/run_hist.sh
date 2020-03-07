#!/bin/bash

for i in 60.5 59.0 57.5 56.0 54.5 53.0 51.5 50.0 48.5 47.0 45.5 44.0 42.5 41.0 39.5 38.0 36.5 35.0 33.5 32.0
do

echo "$i"
./generate_hist.py -i ../${i}/dist.dat > hist_${i}.dat

done
