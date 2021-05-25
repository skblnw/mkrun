#!/bin/bash

filename="result/bm-result"
rm -f $filename
touch $filename
for ii in 2 4 8 16 32
do
    total=0
    for jj in $(seq 1 3); do
      ns=`grep "Performance" result/${ii}-${jj} | awk '{print $2}'`
      total=`echo $total + $ns | bc`
    done
    total=`echo "$total / $jj" | bc -l`
    printf "%2d %6.2f\n" $ii $total >> $filename
done