#!/bin/bash

if [ $# -eq 0 ]; then
    echo "mknamd_prepare <mini|heat|pre|cons|npt1|npt2|smd|us1|us2"
    exit 0
fi

TEMPLATE_DIR=/home/kevin/Dropbox/QWD/scripts/MD/NAMD/2.10
TEMPLATE=$TEMPLATE_DIR/template-make-prepare-sh
cat $TEMPLATE > make_prepare.sh
cp -v $TEMPLATE_DIR/{template-namd,make_submit.sh,vm_getcell.tcl,vm_writepdb.tcl} .

mkdir -p output log

# Choose for cluster gonna use
PS3='Please enter your choice: '
options=("Combo" "Titan" "College" "Let me decide later")
select name in "${options[@]}"
do
    case $name in
        "Combo")
            echo "you chose choice 1"
            sed -i 's/^combo=.*$/combo=true/g' make_prepare.sh
            cp -v $TEMPLATE_DIR/Default/Combo/* .
            break
            ;;
        "Titan")
            echo "you chose choice 2"
            sed -i 's/^titan=.*$/titan=true/g' make_prepare.sh
            cp -v $TEMPLATE_DIR/Default/Titan/* .
            break
            ;;
        "College")
            echo "you chose choice 3"
            sed -i 's/^college=.*$/college=true/g' make_prepare.sh
            break
            ;;
        "Let me decide later")
            break
            ;;
        *) echo invalid option;;
    esac
done

# Ask how many nodes
read -p "How many node you wanna use? `echo $'\n> '`" nnode
if ! [[ $nnode =~ ^-?[0-9]+$ ]]; then
    echo "Must be an integer!"
    exit 1
else
    sed -i 's/^nnode=.*$/nnode='$nnode'/g' make_prepare.sh
fi

# Ask how many cpu per node
read -p "How many cpu are there in 1 node? [16] `echo $'\n> '`" ncpu
if ! [[ $ncpu =~ ^-?[0-9]+$ ]]; then
    echo "Must be an integer!"
    exit 1
else
    sed -i 's/^ncpu=.*$/ncpu='$ncpu'/g' make_prepare.sh
fi

# Turn on options in make_prepare.sh
for name in $@; do
    sed -i 's/^'$name'=.*$/'$name'=true/g' make_prepare.sh
    case $name in
        "us1")
            echo "Info: In case of [$name], we copy some more templates"
            cp -v $TEMPLATE_DIR/US/template-* .
            cp -v $TEMPLATE_DIR/US/*.tcl .
            break
            ;;
        "us2")
            echo "Info: In case of [$name], we copy some more templates"
            cp -v $TEMPLATE_DIR/US/template-* .
            cp -v $TEMPLATE_DIR/US/*.tcl .
            break
            ;;
        *) ;;
    esac
done

# Ask if protein-only
read -p "Is it a protein-only simulation?" -n 1 -r
echo 
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Info: Turn off FlexibleCell"
    sed -i 's/^    useFlexibleCell      yes/    useFlexibleCell      no/g' template-namd
    echo "Info: Turn off ConstantRatio"
    sed -i 's/^    useConstantRatio     yes/    useConstantRatio     no/g' template-namd
fi
