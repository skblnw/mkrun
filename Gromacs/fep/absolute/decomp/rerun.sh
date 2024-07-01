#!/bin/bash

for ii in $(seq 0 47); do
    DIR=/data/kevin/1uw6_chargedlig6/Lambda_$ii/Production_MD
    jj=$((ii+1))

    sed -i "s/init_lambda_state.*/init_lambda_state = $ii/g" win0.mdp
    gmx grompp -o $DIR/win$ii.tpr -f win0.mdp -c $DIR/md.gro -r $DIR/md.gro -p topol.top -n index.ndx -t $DIR/md.cpt -maxwarn 1
    gmx mdrun -ntmpi 1 -ntomp 16  -s $DIR/win$ii.tpr -rerun $DIR/md.trr -deffnm $DIR/win$ii

    sed -i "s/init_lambda_state.*/init_lambda_state = $jj/g" win1.mdp
    gmx grompp -o $DIR/win$jj.tpr -f win1.mdp -c $DIR/md.gro -r $DIR/md.gro -p topol.top -n index.ndx -t $DIR/md.cpt -maxwarn 1
    gmx mdrun -ntmpi 1 -ntomp 16 -s $DIR/win$jj.tpr -rerun $DIR/md.trr -deffnm $DIR/win$jj
done
