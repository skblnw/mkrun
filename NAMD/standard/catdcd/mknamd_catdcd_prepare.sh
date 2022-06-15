#!/bin/bash
#########################################
## Description: Utilize CatDCD binary in a smart way.
## 1. By default, it cat all dcd (1 to N) you input into one single dcd
## 2. It would check if you do not start from 1, whether an old one with index-1 exists
## 3. One new tmp.dcd will be cat and app to the old one
## Author: Kev, Jul 2016
## Usage: make_XXX </path/to/trajectory/NPT-{1..N}.dcd>
## Input:
## Output:
## Units: 
## Other Notes: 
#########################################

if [ "$#" -eq 0 ] || [ "$#" -gt 1 ]; then
    echo "make_XXX </path/to/trajectory/{mini/heat/pre/cons}.dcd>"
    exit 1
fi

SRC_DIR=$1
OUTPUT_DIR=`pwd`

read -p "[Old DCD] How many ps per frame? [1000] " ps
if ! [[ $ps =~ ^-?[0-9]+$ ]]; then
    echo "Must be an integer!"
    exit 1
fi

read -p "[Input DCD] How many ps per frame? [2fs*1000=2] " ps2
if ! [[ $ps2 =~ ^-?[0-9]+$ ]]; then
    echo "Must be an integer!"
    exit 1
fi

read -p "Name of the index-XXX.txt? [noW] " ind

stride=$(( $ps / $ps2 ))
echo -e "The stride is $stride"
cd $SRC_DIR
catdcd -o $OUTPUT_DIR/preparation-pf${ps}ps.dcd -i $OUTPUT_DIR/index-$ind.txt -stride $stride mini-*.dcd heat.dcd
