#!/bin/bash

previous=step4_eq_npt
prefix=step5_md
for NN in {61..70}; do
    MM=$(expr $NN - 1)
    if [ $MM -eq 0 ]; then
        sed -e 's/^iname=.*$/iname='$previous'/g' \
            -e 's/^oname=.*$/oname='$prefix'-'$NN'/g' \
            -e 's/PREFIX/'$prefix'/g' \
            -e 's/NN/'$NN'/g' \
            -e 's/^CONTINUE=.*$/CONTINUE=false/g' \
            template-job > job_$prefix-$NN.pbs
    else
        sed -e 's/^iname=.*$/iname='$prefix'-'$MM'/g' \
            -e 's/^oname=.*$/oname='$prefix'-'$NN'/g' \
            -e 's/PREFIX/'$prefix'/g' \
            -e 's/NN/'$NN'/g' \
            -e 's/^CONTINUE=.*$/CONTINUE=true/g' \
            template-job > job_$prefix-$NN.pbs
    fi
done
