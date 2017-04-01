#!/bin/bash

# Some Global variables
TEMPLATE_DIR=./templates

cp -f templates/template-mart-pbs job-mart.pbs

while true; do

    read -p "How many cores?" REPLY </dev/tty 
    
    if [[ "$REPLY" =~ ^[0-9]+$ ]]; then
        sed -i 's/NPPN=.*$/NPPN='$REPLY'/g' job-mart.pbs
        break
    fi

done

if ask "mini?"; then
    sed -i 's/^mini=.*$/mini=true/g' job-mart.pbs

    if ask "Choose steep as minimization scheme?`echo $'\n '`[y:Steep|n:CG]"; then
        cp templates/template-mart-mini-steep-mdp mini.mdp
    else
        cp templates/template-mart-mini-cg-mdp mini.mdp
    fi
fi

if ask "eq?"; then
    sed -i 's/^eq=.*$/eq=true/g' job-mart.pbs

    cp templates/template-mart-eq-npt-mdp eq.mdp
fi

if ask "md1?"; then
    sed -i 's/^md1=.*$/md1=true/g' job-mart.pbs

    cp templates/template-mart-md-npt-mdp md.mdp
fi

if ask "md2?"; then
    sed -i 's/^md2=.*$/md2=true/g' job-mart.pbs
fi
