#!/bin/bash

for ii in {11..11}
do
    jj=$(expr $ii - 1)
    sed -e 's/MM/'$jj'/g' -e 's/NN/'$ii'/g' template-namd > NPT-${ii}.namd
    sed -e 's/MM/'$jj'/g' -e 's/NN/'$ii'/g' template-pbs > pitzer.pbs
    LAST=$(qsub -W depend=afterany:21225 pitzer.pbs)
done

for ii in {12..20}
do
    jj=$(expr $ii - 1)
    sed -e 's/MM/'$jj'/g' -e 's/NN/'$ii'/g' template-namd > NPT-${ii}.namd
    sed -e 's/MM/'$jj'/g' -e 's/NN/'$ii'/g' template-pbs > pitzer.pbs
    LAST=$(qsub -W depend=afterany:$LAST pitzer.pbs)
done
