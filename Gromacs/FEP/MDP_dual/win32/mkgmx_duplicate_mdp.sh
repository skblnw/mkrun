#!/bin/bash

if [[ "$#" -ne 1 ]]; then
    echo "Usage: mkgmx_duplicate_mdp input.mdp"
    exit 0
fi

TEMPLATE="$1"

prefix=$(basename $TEMPLATE .mdp)
echo "mkgmx> The prefix is $prefix"

for ii in $(seq 0 34)
do
    echo "mkgmx> Initial Lambda State $ii"

    sed -e 's/^init_lambda_state.*$/init_lambda_state = '$ii'/g' \
    $TEMPLATE > ${prefix}_${ii}.mdp
done
