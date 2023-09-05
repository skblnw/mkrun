#!/bin/bash

jj=0
for ii in 50.0 48.5 47.0 45.5
do
    export CUDA_VISIBLE_DEVICES=$jj

    cp *.in COM_dist.RST $ii
    cd $ii
    sed -i "s/DISTHERE/${ii}/g" COM_dist.RST
    bash ../run.sh &
    cd ..

    jj=$((jj+1))
done
