#!/bin/bash

jj=0
for ii in 44.0 42.5 41.0 39.5
do
    export CUDA_VISIBLE_DEVICES=$jj

    cp *.in $ii
    cd $ii
    bash ../run.sh &
    cd ..

    jj=$((jj+1))
done
