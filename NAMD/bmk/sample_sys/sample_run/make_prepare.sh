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
# 17 May 2016
# -Make it a template for preparing any MD program
################################################
# NAMD Commands (GPU):
# COMBO: \$CHARMRUN +p\$NP2 ++ppn 15 ++scalable-start ++nodelist \$nodefile \$NAMD +setcpuaffinity +pemap 1-15 +commap 0 ++verbose
# Titan: aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0
# College: mpirun -n \$NP \$NAMD
################################################
# NAMD Commands (CPU):
# COMBO: \$CHARMRUN +p\$NP ++ppn 16 ++scalable-start ++nodelist \$nodefile \$NAMD ++verbose
################################################
# Gromacs Commands (GPU):
# COMBO: $MPIRUN -mca btl self,openib -np $NN -npernode 1 $GMXBIN mdrun -v -ntomp 7 -pin on -s W50k.tpr -deffnm output/gpu-1node-NNomp
################################################
#
# Steps:
# 1. copy command above into variable COMMAND
# 2. change $CHARMRUN and $NAMD to $NAMD_PATH/charmrun
# 3. change nnode list
# 4. EXECUTE!


nheader=24

bmk=true
bmk2=false

if $bmk; then
    cpu=15
    for nnode in {1,2,3,4,8,12}
    do
        ncpu=$(($cpu * $nnode))
        NAMD_PATH="/share/apps/NAMD_2.11_Linux-x86_64-ib-smp-cuda"
        COMMAND="$NAMD_PATH/charmrun +p$ncpu ++ppn 15 ++nodelist \$nodefile $NAMD_PATH/namd2 +setcpuaffinity +pemap 1-15 +commap 0"
        #echo "#!/bin/bash" > job-1-$ncpu.pbs
        sed -e 's/^#PBS -N .*/#PBS -N bmk-'$nnode'/g' \
            -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn=16/g' \
            template-job-pbs > job-$nnode.pbs
        sed -i $nheader',$d' job-$nnode.pbs

        sed -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME bmk-'$nnode'/g' \
            -e 's/^set TS.*$/set TS 1000/g' \
            template-namd > outputNN/bmk-$nnode.namd
        echo -e "$COMMAND outputNN/bmk-$nnode.namd > NN/bmk-$nnode.log\nwait" >> job-$nnode.pbs
    done
fi


if $bmk2; then
    for ncpu in {2,4,6,8,12,16}
    do
        nppn=$(($ncpu - 1))
        NAMD_PATH="/share/apps/NAMD_2.11_Linux-x86_64-ib-smp-cuda"
        COMMAND="$NAMD_PATH/charmrun +p\$NP2 ++ppn $nppn ++nodelist \$nodefile $NAMD_PATH/namd2 +setcpuaffinity +pemap 1-$nppn +commap 0"
        name=12-$ncpu
        #echo "#!/bin/bash" > job-1-$ncpu.pbs
        sed -e 's/^#PBS -N .*/#PBS -N bmk-'$name'/g' \
            -e 's/^#PBS -l nodes=.*/#PBS -l nodes=12:ppn='$ncpu'/g' \
            template-job-pbs > job-$name.pbs
        sed -i $nheader',$d' job-$name.pbs

        sed -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME bmk-12-'$ncpu'/g' \
            -e 's/^set TS.*$/set TS 2000/g' \
            template-namd > outputNN/bmk-$name.namd
        echo -e "$COMMAND outputNN/bmk-$name.namd > NN/bmk-$name.log\nwait" >> job-$name.pbs
    done
fi
