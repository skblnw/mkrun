#!/bin/bash
#########################################
## Description: This NAMD Benchmark Grepping script is for daily usage of 
## 1. finding out all "# of cores used" on the Benchmark line of an NAMD log output
## 2. sort and average over same "#"
## Author: Kev, Jul 2016
## Usage: 
## Input:
## Output:
## Units: 
## Other Notes: 
#########################################

if [ $# -eq 0 ]; then
    echo "make_bmk-routine <*.log>"
    exit 1
fi

ncpu=`grep "Benchmark time" $@ | awk '{print $4}' | sort -u -h | xargs`
echo -e "# of cores found: $ncpu \nNow compute benchmarks"

filename="bm-result"
if [ -s $filename.dat ]; then
    read -p "$filename.dat already exists! Overwite? [y/N]" -n 1 -r
    echo 
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    else
        rm -f $filename.dat
    fi
fi

echo -e "CPU#\tNode#\tBM\tN" > $filename.dat

for ii in $ncpu; do
    DAY=`grep "Benchmark time" $@ | awk '{print $4" "$8}' | grep "^$ii " | awk '{sum+=$2} END {print sum/NR}'`
    NS=`printf '%.2f\n' "$(echo "scale=2; 1.0/$DAY" | bc)"`
    NLINE=`grep "Benchmark time" $@ | awk '{print $4" "$8}' | grep "^$ii " | wc -l`
    NLOG=`echo "$NLINE/6" | bc`
    NODE=`grep 'physical nodes' $@ | grep " $ii processors" | awk '{print $8}' | uniq`
    echo -e "$ii\t$NODE\t$NS\t$NLOG" >> $filename.dat
done
