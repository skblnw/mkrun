#!/bin/bash

for ii in 1; do
filename=smd$ii
LOG=/share/data/kevin/PH-memb/smd/4/run-ss-xy/log/$filename-1.log
RMSD_PSPF=10
SEP_NF=$(($RMSD_PSPF / 2))

grep "^SMD " $LOG | awk -v nth="$SEP_NF" 'NR % nth == 0 {print $5}' > namdlog/$filename.dat
done
