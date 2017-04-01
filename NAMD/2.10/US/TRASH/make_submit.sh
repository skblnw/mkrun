LAST=$(qsub job-1.pbs)
echo $LAST

for ii in {2..5}
do
    LAST=$(qsub -W depend=afterany:$LAST job-$ii.pbs)
    echo $LAST
done
