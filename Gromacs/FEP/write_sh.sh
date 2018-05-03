#!/bin/bash

prefix=9state
last=md\$LAMBDA
next=md\$LAMBDA-2

for ii in $(seq 0 8)
do
    sed -e 's/^LAMBDA=.*$/LAMBDA='$ii'/g' \
    -e 's/^LAST=.*$/LAST='$last'/g' \
    -e 's/^NEXT=.*$/NEXT='$next'/g' \
    template_continuation > job-${prefix}_${ii}.sh
done