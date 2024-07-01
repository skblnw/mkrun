#!/bin/bash

prefix=20state

dir_mdp=MDP_20

for ii in $(seq 0 20)
do
    echo ">  Initial Lambda State $ii"

    sed -e 's/^init_lambda_state.*$/init_lambda_state = '$ii'/g' \
    $dir_mdp/template_md > $dir_mdp/md$ii.mdp

    sed -e 's/^LAMBDA=.*$/LAMBDA='$ii'/g' \
    template_job > job-${prefix}_${ii}.sh
done


# Continuation
if false; then
    
last=md\$LAMBDA
next=md\$LAMBDA-2

for ii in $(seq 0 8)
do
    echo ">  Initial Lambda State $ii"

    sed -e 's/^LAMBDA=.*$/LAMBDA='$ii'/g' \
    -e 's/^LAST=.*$/LAST='$last'/g' \
    -e 's/^NEXT=.*$/NEXT='$next'/g' \
    template_continuation > job-${prefix}_${ii}.sh
done

fi