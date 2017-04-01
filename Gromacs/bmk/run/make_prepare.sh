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
# NAMD Commands:
# COMBO: \$CHARMRUN +p\$NP2 ++ppn 15 ++scalable-start ++nodelist \$nodefile \$NAMD +setcpuaffinity +pemap 1-15 +commap 0 ++verbose
# Titan: aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0
# College: mpirun -n \$NP \$NAMD
################################################
# Gromacs Commands:
# Combo: $MPIRUN -mca btl self,openib -np $NN -npernode 1 $GMXBIN mdrun -v -ntomp 7 -pin on -s W50k.tpr -deffnm output/gpu-1node-NNomp
################################################

nheader=26

bmk=false      # Ordinary benchmark
ppn=false      # Process/Core per node test usually for CUDA version of the MD program

MPIRUN=""
GMXBIN=""

if $bmk; then
    ncpu=16
    for nnode in {2,4,6,8}
    do
        name=$nnode
        sed -e 's/^#PBS -N .*/#PBS -N bmk-'$name'/g' \
            -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn=16/g' \
            template-job-pbs > job-$name.pbs
        sed -i $nheader',$d' job-$name.pbs
        np=$(($nnode * 2))
        COMMAND="\$MPIRUN -mca btl self,openib -np $nnode -npernode 1 --hostfile ./nodefile \$GMXBIN mdrun -v -pin on"
        echo -e "$COMMAND -ntomp 8 -s tpr/W50k.tpr -deffnm NUM_ITER/50k-gpu-${nnode}node-8omp\nwait" >> job-$name.pbs
    done
fi

if $ppn; then
    for ncpu in {1,2,4,8,12,16}
    do
        nnode=4
        name=$ncpu
        sed -e 's/^#PBS -N .*/#PBS -N bmk-'$name'/g' \
            -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn=16/g' \
            template-job-pbs > job-$name.pbs
        sed -i $nheader',$d' job-$name.pbs
        np=$(($nnode * 2))
        COMMAND="\$MPIRUN -mca btl self,openib -np $nnode -npernode 1 --hostfile ./nodefile \$GMXBIN mdrun -v -pin on"
		#COMMAND="\$GMXBIN mdrun -v -pin on"
		# Usually MD programs on one single node behave very differetly compared to MPI. So test explicitly.
        echo -e "$COMMAND -ntomp $ncpu -s tpr/W50k.tpr -deffnm NUM_ITER/50k-gpu-${nnode}node-${ncpu}omp\nwait" >> job-$name.pbs
    done
fi