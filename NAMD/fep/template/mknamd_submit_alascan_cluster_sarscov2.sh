#!/bin/bash

START=2
REPLICA=5
RESID="339 371 373 375 417 440 446 477 478 493 496 498 501 505"

[ $# -eq 0 ] && { echo "mkpbs> Usage: $0 [NJOB] [-prepare|-submit]"; echo "mkpbs> NJOB: Number of jobs for each mutation (Default: 0). If 0, then $REPLICA runs will be executed inside each job."; exit 1; }
OPT="$2"

for id in $RESID; do
    for prefix in bound; do
        cd posi$id/$prefix
        ln -sfn ~/toppar_c36_jul20/ toppar
        ln -sfn ~/toppar_water_ions_namd.str .
        if [ $1 -gt 0 ]; then
            for jj in $(seq 1 $1); do
                cat > pbs$jj <<EOF
#PBS -N $id-$prefix
#PBS -q fep
#PBS -l nodes=1:ppn=1:gpus=1
#PBS -l walltime=12:00:00
#PBS -S /bin/bash
#PBS -j oe
date
hostname
echo "CUDA_VISIBLE_DEVICES: \$CUDA_VISIBLE_DEVICES"
export LD_LIBRARY_PATH=/public/software/lib/:\$LD_LIBRARY_PATH
source /public/software/compiler/intel/intel-compiler-2017.5.239/bin/compilervars.sh intel64
cd \$PBS_O_WORKDIR
echo \$PBS_O_WORKDIR
NAMD="/public/software/apps/NAMD_3.0alpha9/namd3"
\$NAMD +p4 +devices 0 eq/fep.eq.namd >& eq/LOG_eq
rsync -a eq/fep.namd eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc trial$jj
\$NAMD +p1 +devices 0 trial$jj/fep.namd >& trial$jj/LOG_fep
EOF
                case $OPT in
                    -prepare)
                    echo "mkpbs> Created PBS for [Position $id] [$prefix state]"
                    ;;
                    -submit)
                    echo "mkpbs> Submitted for [Position $id] [$prefix state]"
                    qsub pbs$jj
                    ;;
                    *)    # unknown option
                    echo -e "Unknown option $OPT, must be either -prepare or -submit"; exit 0
                    ;;
                esac
            done
        else
                cat > pbs <<EOF
#PBS -N $id-$prefix
#PBS -q fep
#PBS -l nodes=1:ppn=1:gpus=1
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -j oe
date
hostname
echo "CUDA_VISIBLE_DEVICES: \$CUDA_VISIBLE_DEVICES"
export LD_LIBRARY_PATH=/public/software/lib/:\$LD_LIBRARY_PATH
source /public/software/compiler/intel/intel-compiler-2017.5.239/bin/compilervars.sh intel64
cd \$PBS_O_WORKDIR
echo \$PBS_O_WORKDIR
NAMD="/public/software/apps/NAMD_3.0alpha9/namd3"
for ii in \$(seq $START $REPLICA); do
    rsync -a eq/fep.namd eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc trial\$ii
    \$NAMD +p1 +devices 0 trial\$ii/fep.namd >& trial\$ii/LOG_fep
done
date
EOF
            case $OPT in
                -prepare)
                echo "mkpbs> Created PBS for [Position $id] [$prefix state]"
                ;;
                -submit)
                echo "mkpbs> Submitted for [Position $id] [$prefix state]"
                qsub pbs
                ;;
                *)    # unknown option
                echo -e "Unknown option $OPT, must be either -prepare or -submit"; exit 0
                ;;
            esac
        fi
        cd ../..
    done
done
