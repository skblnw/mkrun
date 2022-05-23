#!/bin/bash

[ $# -ne 2 ] && { echo "mknamd> Usage: $0 [start] [end]"; exit 1; }

START=$1
END=$2
RUN_DIR=run
NORMRUN="namd3 +p2 +devices 0"
CUDARUN="namd3 +p1 +devices 0"

add_pbs () {
    if [ "$#" -ne 2 ]; then
        echo "> Wrong argument number"
        exit 1
    fi

    echo "echo -e \"mknamd> Working hard with $1...\nmknamd> Finished at:\"" >> $2
    if $acceleration; then
        echo -e "$CUDARUN $1.namd > log/$1.log\nwait\ndate" >> $2
    else
        echo -e "$NORMRUN $1.namd > log/$1.log\nwait\ndate" >> $2
    fi
}

rm -r $RUN_DIR/run.sh
acceleration=false
for ii in $(seq $START $END); do
    prefix=md
    frequency=1000

    if [ $ii -eq 1 ]; then
        inputname="output\/cons-3"
        outputname="output\/${prefix}.${ii}"
    else
        jj=$((ii-1))
        inputname="output\/${prefix}.${jj}"
        outputname="output\/${prefix}.${ii}"
    fi

    sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
        -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.restart.vel/g' \
        -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.restart.coor/g' \
        -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.restart.xsc/g' \
        -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
        -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
        -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
        -e 's/^COMmotion.*$/COMmotion no/g' \
        -e 's/^set ITEMP.*$/set ITEMP 303/g' \
        -e 's/^set FTEMP.*$/set FTEMP 303/g' \
        -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
        -e 's/^set TS.*$/set TS 500000/g' \
        -e 's/^set CUDASOA.*$/set CUDASOA 0/g' \
        -e 's/^set CONSSCALE.*$/set CONSSCALE 10/g' \
        -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/cons_posres/g' \
        -e 's/^set COLFILE.*$/set COLFILE ..\/umbrella.col/g' \
        template-namd > $RUN_DIR/${prefix}.${ii}.namd

    add_pbs ${prefix}.${ii} $RUN_DIR/run.sh
done
