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

if [ $# -eq 0 ]; then
    echo "make_XXX </path/to/trajectory/NPT-{1..N}.dcd>"
    exit 1
fi

SKIP_OLD=0

prefix=`basename "${@: -1}" | awk -F '-' '{print $1}'`
prefix_lower=`echo "$prefix" | tr '[:upper:]' '[:lower:]'`

FIRST=`echo $1 | sed -n -e 's/.*'$prefix'-//p' | sed -n -e 's/.dcd//p'`
if ! [[ $FIRST =~ ^-?[0-9]+$ ]]; then
    echo "'FIRST' Must be an integer!"
    read -p "Type your 'FIRST' " FIRST
    if ! [[ $FIRST =~ ^-?[0-9]+$ ]]; then
        echo "Must be an integer!"
        exit 1
    fi
    SKIP_OLD=1
fi

LAST=`echo "${@: -1}" | sed -n -e 's/.*'$prefix'-//p' | sed -n -e 's/.dcd//p'`
if ! [[ $LAST =~ ^-?[0-9]+$ ]]; then
    echo "'LAST' Must be an integer!"
    read -p "Type your 'LAST' " LAST
    if ! [[ $LAST =~ ^-?[0-9]+$ ]]; then
        echo "Must be an integer!"
        exit 1
    fi
fi

if [ "$FIRST" -gt 2 ] && [ "$SKIP_OLD" -eq 0 ]; then
    echo -e "Smallest DCD index larger than 2\nTry to find an old dcd"

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

    LLAST=$(( $FIRST - 1))
    read -p "Name of the index-XXX.txt? [noW] " ind
    filename=$prefix_lower-$ind-pf${ps}ps
    if [ ! -s $filename-$LLAST.dcd ]; then
        echo -e "$filename-$LLAST.dcd does not exist!"
        exit 1
    else
        echo -e "Found "$filename-$LLAST.dcd"!"
    fi

    if [ -s $filename-$LAST.dcd ]; then
        echo -e "$filename-$LAST.dcd already exist!"
        exit 1
    fi

    stride=$(( $ps / $ps2 ))
    echo -e "The stride is $stride"
    catdcd -o $filename-new.dcd -i index-$ind.txt -stride $stride $@
    catdcd -o $filename-$LAST.dcd $filename-$LLAST.dcd $filename-new.dcd
    rm -f $filename-new.dcd
else
    echo -e "DCD index starts from $FIRST. Make you a trimmed one :)"
    read -p "[Input DCD] How many ps per frame? [2fs*1000=2] " ps2
    if ! [[ $ps2 =~ ^-?[0-9]+$ ]]; then
        echo "Must be an integer!"
        exit 1
    fi
    read -p "[Output DCD] How many ps per frame? [1000] " ps
    if ! [[ $ps =~ ^-?[0-9]+$ ]]; then
        echo "Must be an integer!"
        exit 1
    fi
    read -p "Name of the index-XXX.txt? [noW] " ind
    filename=$prefix_lower-$ind-pf${ps}ps
    stride=$(( $ps / $ps2 ))
    echo -e "The stride is $stride"
    catdcd -o $filename-$LAST.dcd -i index-$ind.txt -stride $stride $@
fi
