#!/bin/bash
#########################################
## Description: Make script for BMW water model MD parameters
## Author: Kevin@CityUHK Mar 2016
## Usage: 
## Input:
## Output:
## Units: 
## Other Notes: Modify from make script for Martini MD parameters
#########################################

# Some Global variables
TEMPLATE_DIR=./templates # folder containing all templates

cp -f $TEMPLATE_DIR/template-bmw-pbs job-bmw.pbs
cp -f $TEMPLATE_DIR/table*.xvg .


# Ask for number of nodes or cores
while true; do
    read -p "How many <nodes>? [integer|Press Enter] " NN </dev/tty 
    if [[ "$NN" =~ ^[0-9]+$ ]]; then
        read -p "How many <cores> per node? " NPPN </dev/tty 
        if [[ "$NPPN" =~ ^[0-9]+$ ]]; then
            sed -i 's/^#PBS -l nodes=1:ppn=16/#PBS -l nodes='$NN':ppn='$NPPN'/g' job-bmw.pbs
            NP=$(($NPPN * $NN))
            sed -i 's/NP=.*$/NP='$NP'/g' job-bmw.pbs
            break
        fi
    else
        read -p "So you wanna specify # of <cores> in total? How many? " REPLY </dev/tty 
        if [[ "$REPLY" =~ ^[0-9]+$ ]]; then
            sed -i 's/NP=.*$/NP='$REPLY'/g' job-bmw.pbs
            break
        fi
    fi
done

if ask "flexible?" N; then
    if [[ -s ../system.top ]]; then
        sed -e 's/martini_v2.1_bmw.itp/martini_v2.1_flexible.itp/g' ../system.top > ../flexible.top
	else
		echo "You don't have the system.top file!"
		exit 0
	fi

    sed -i 's/^flexible=.*$/flexible=true/g' job-bmw.pbs

    cp -f $TEMPLATE_DIR/template-bmw-mini-mdp bmw_mini.mdp
fi

if ask "mini?" N; then
    sed -i 's/^mini=.*$/mini=true/g' job-bmw.pbs

    cp -f $TEMPLATE_DIR/template-bmw-mini-mdp bmw_mini.mdp
fi

if ask "nvt?" N; then
    sed -i 's/^nvt=.*$/nvt=true/g' job-bmw.pbs

    cp -f $TEMPLATE_DIR/template-bmw-nvt-mdp bmw_nvt.mdp
fi

if ask "npt?" N; then
    sed -i 's/^npt=.*$/npt=true/g' job-bmw.pbs

    cp -f $TEMPLATE_DIR/template-bmw-npt-mdp bmw_npt.mdp
fi

if ask "md1?" N; then
    sed -i 's/^md1=.*$/md1=true/g' job-bmw.pbs

    cp -f $TEMPLATE_DIR/template-bmw-npt-mdp bmw_md.mdp
fi

if ask "md2?" N; then
    sed -i 's/^md2=.*$/md2=true/g' job-bmw.pbs

    while true; do
        read -p "How many <ps> do you wanna extend? " NSTEP </dev/tty
        if [[ "$NSTEP" =~ ^[0-9]+$ ]]; then
            sed -i 's/-extend.*0/-extend '$NSTEP'/g' job-bmw.pbs
            break
        fi
    done

    while true; do
        read -p "Starting index? " STARTJ </dev/tty
        if [[ "$STARTJ" =~ ^[0-9]+$ ]]; then
            read -p "Ending index? " ENDJ </dev/tty
            if [[ "$ENDJ" =~ ^[0-9]+$ ]]; then
                sed -i -e 's/STARTJ/'$STARTJ'/g' -e 's/ENDJ/'$ENDJ'/g' job-bmw.pbs
                break
            fi
        fi
    done
fi
