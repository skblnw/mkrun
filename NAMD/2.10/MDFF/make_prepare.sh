################################################
# 11 Sep 2015
# CBMSG
# Kevin
# I have prepared this script solely for preparing NAMD MDFF job scripts
#
################################################
# MD Commands:
# COMBO: \$CHARMRUN +p\$NP2 ++ppn 15 ++scalable-start ++nodelist \$nodefile \$NAMD +setcpuaffinity +pemap 1-15 +commap 0 ++verbose
# Titan: aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0
# College: mpirun -n \$NP \$NAMD
################################################

nnode=60
ncpu=15
nheader=33
COMMAND="aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0"

heat=true
run1=true
run2=true

if $heat; then
    jobname=heat
    sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
        -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
        -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
        template-job-pbs > job-$jobname.pbs
    sed -i $nheader',$d' job-$jobname.pbs

    sed -e 's/^set INPUTNAME.*$//g' \
	-e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/mdff-heat/g' \
        -e 's/^set GSCALE.*$/set GSCALE 0.1/g' \
	-e 's/^set ITEMP.*$/set ITEMP 0/g' \
        -e 's/^set MS.*$/set MS 5000/g' \
        -e 's/^set TS.*$/set TS 1000000/g' \
	-e 's/^symmetryRestraints.*$/symmetryRestraints off/g' \
        template-mdff-namd > mdff-heat.namd
    echo -e "$COMMAND mdff-heat.namd > mdff-heat.log\nwait" >> job-$jobname.pbs
fi

if $run1; then
    jobname=run1
    sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
        -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
        -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
        template-job-pbs > job-1.pbs
    sed -i $nheader',$d' job-1.pbs

    previous=heat
    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/mdff-'$previous'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/mdff-1/g' \
        -e 's/^set GSCALE.*$/set GSCALE 0.1/g' \
	-e 's/^set MS.*$/set MS 0/g' \
        -e 's/^set TS.*$/set TS 1000000/g' \
	-e 's/^symmetryfirststep.*$/symmetryfirststep 1/g' \
	-e 's/^symmetryfirstfullstep.*$/symmetryfirstfullstep 1000000/g' \
        template-mdff-namd > mdff-1.namd
    echo -e "$COMMAND mdff-1.namd > mdff-1.log\nwait" >> job-1.pbs
fi

if $run2; then
    for ii in {2..2}
    do
        jobname=run$ii
        sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
            -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
            -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
            template-job-pbs > job-$ii.pbs
        sed -i $nheader',$d' job-$ii.pbs

        jj=$(expr $ii - 1)
        sed \
            -e 's/^set INPUTNAME.*$/set INPUTNAME output\/mdff-'$jj'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/mdff-'$ii'/g' \
            -e 's/^set GSCALE.*$/set GSCALE 0.1/g' \
            -e 's/^set MS.*$/set MS 0/g' \
            -e 's/^set TS.*$/set TS 1000000/g' \
            -e 's/^symmetryfirststep.*$/symmetryfirststep 1/g' \
            -e 's/^symmetryfirstfullstep.*$//g' \
            template-mdff-namd > mdff-$ii.namd
        echo -e "$COMMAND mdff-$ii.namd > mdff-$ii.log\nwait" >> job-$ii.pbs
    done
fi
