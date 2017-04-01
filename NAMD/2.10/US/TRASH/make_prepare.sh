################################################
# 4 Aug 2015
# CBMSG
# Kevin
# I have prepared this script solely for preparing NAMD MD job scripts universally
#
# Updates:
# 13 Sep 2015
# -Added option bmk
# 4 Aug 2015
# -Added variable previous
# 3 Aug 2015
# -Added variable COMMAND
# -Added variable jobname
# -Added variable nnode and nheader
################################################
# MD Commands:
# COMBO: \$CHARMRUN +p\$NP2 ++ppn 15 ++scalable-start ++nodelist \$nodefile \$NAMD +setcpuaffinity +pemap 1-15 +commap 0 ++verbose
# Titan: aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0
# College: mpirun -n \$NP \$NAMD
################################################

tnode=264
nnode=24
ncpu=16
nheader=36
COMMAND="aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0"

npt2=true

if $npt2; then
    for ii in {1..1}
    do
        jobname=us-$ii
        sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
            -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$tnode'/g' \
            -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=12:00:00/g' \
            template-job-pbs-titan > job-$ii.pbs
        sed -i $nheader',$d' job-$ii.pbs

        for dd in 860 875 890 905 920 935 950 965 980 990 1000
        do
            aa=`echo "scale=1;$dd / 10" | bc`  # compute actual position of window centers in Angstron
            jj=$(expr $ii - 1)
            sed -e 's/NN/'$ii'/g' -e 's/MM/'$jj'/g' -e 's/DD/'$dd'/g' template-us-namd > us-z${dd}-${ii}.namd
            sed -e 's/POSaa/'$aa'/g' template-colvars > colvars-z$dd.in
#            sed -e 's/POSaa/'$dd'/g' -e 's/NUM/'$ii'/g' template-gjob > job-z${dd}-${ii}.pbs
#            sed -e 's/POSaa/'$dd'/g' -e 's/NUM/'$ii'/g' template-job-pbs-combo > job-z${dd}-${ii}.pbs
            echo -e "$COMMAND us-z${dd}-${ii}.namd > us-z${dd}-${ii}.log &" >> job-$ii.pbs
        done

        echo -e "wait" >> job-${ii}.pbs
    done
fi
