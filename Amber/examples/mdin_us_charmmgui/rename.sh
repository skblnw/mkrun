#!/bin/bash

for ii in $(seq 1 32)
do
    mkdir $ii
    mv ionized.parm7.$ii $ii
    mv ionized.rst7.$ii $ii
done
