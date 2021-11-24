#!/bin/bash

CUDARUN="namd3 +p1 +devices 0"
md_posres=false

check_exist () {
    if [ ! -s "$1" ]; then
        echo "> $1 DOES NOT EXIST! Be careful!"
        exit 1
    fi
}

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

previous="md"
prefix="t1"

check_exist template-namd
check_exist run/output/${previous}.restart.coor
check_exist run/output/${previous}.restart.vel
check_exist run/output/${previous}.restart.xsc
if $md_posres; then
    check_exist restraints/cons_posres.pdb
fi
frequency=50000
acceleration=true
mkdir -p run/log
[ -s "run/lrun.sh" ] && { rm run/lrun.sh; }
for ii in $(seq 2 3)
do

    if [ $ii -eq 2 ]; then
        inputname="output\/${previous}"
        outputname="output\/${prefix}.part${ii}"
    else
        jj=$((ii-1))
        inputname="output\/${prefix}.part${jj}"
        outputname="output\/${prefix}.part${ii}"
    fi

    if $md_posres; then
        sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
            -e 's/^set CONSSCALE.*$/set CONSSCALE 10/g' \
            -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/cons_posres/g' \
            -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.restart.vel/g' \
            -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.restart.coor/g' \
            -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.restart.xsc/g' \
            -e 's/^COMmotion.*$/COMmotion yes/g' \
            -e 's/^set ITEMP.*$/set ITEMP 303/g' \
            -e 's/^set FTEMP.*$/set FTEMP 303/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
            -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
            -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
            -e 's/^set TS.*$/set TS 250000000/g' \
            -e 's/^set CUDASOA.*$/set CUDASOA 1/g' \
            template-namd > run/${prefix}.part${ii}.namd
    else
        sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
            -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.restart.vel/g' \
            -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.restart.coor/g' \
            -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.restart.xsc/g' \
            -e 's/^COMmotion.*$/COMmotion yes/g' \
            -e 's/^set ITEMP.*$/set ITEMP 303/g' \
            -e 's/^set FTEMP.*$/set FTEMP 303/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
            -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
            -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
            -e 's/^set TS.*$/set TS 250000000/g' \
            -e 's/^set CUDASOA.*$/set CUDASOA 1/g' \
            template-namd > run/${prefix}.part${ii}.namd
    fi
    
    add_pbs ${prefix}.part${ii} run/lrun.sh

done
