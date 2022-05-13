#!/bin/bash

NSTEP=350000000
PREFIX=2

[ $# -eq 0 ] && { echo -e "Usage: fixpbs [prefix] [-run|-lrun] [run#]"; exit 1; }

PREFIX="$3"

check_exist () {
    if [ ! -s "$1" ]; then
        echo "> $1 DOES NOT EXIST! Be careful!"
        exit 1
    fi
}

case $2 in
    -run)
    echo "mkpbs> Fix PBS for eq run"
    ;;
    -lrun)
    echo "mkpbs> Fix PBS for a $((NSTEP/500000)) ns run"
    ;;
    *)    # unknown option
    echo -e "Unknown option $2, must be either -run or -lrun"; exit 0
    ;;
esac

cat > mdnamd.pbs <<EOF
#PBS -N $1
#PBS -q md
#PBS -l nodes=1:ppn=1:gpus=1
#PBS -l walltime=24:00:00
#PBS -S /usr/bin/csh
#PBS -j oe

date
echo \$SHELL
echo "CUDA_VISIBLE_DEVICES: \$CUDA_VISIBLE_DEVICES"
cd \$PBS_O_WORKDIR
echo \$PBS_O_WORKDIR

setenv LD_LIBRARY_PATH /public/software/lib/
source /public/software/compiler/intel/intel-compiler-2017.5.239/bin/compilervars.csh intel64
set NAMD="/public/software/apps/NAMD_3.0alpha9/namd3"
EOF

check_exist README

if [[ "$2" == "-lrun" ]]; then
  sed -e "s/stepspercycle           20.*//g" \
      -e "s/pairlistsPerCycle       2.*//g" \
      -e "s/^# run/CUDASOAintegrate on/g" \
      -e "s/binCoordinates          \$inputname.coor/binCoordinates \$inputname.restart.coor/g" \
      -e "s/binVelocities           \$inputname.vel/binVelocities \$inputname.restart.vel/g" \
      -e "s/extendedSystem          \$inputname.xsc/extendedSystem \$inputname.restart.xsc/g" \
      -e "s/^numsteps.*/numsteps $NSTEP/g" \
      -e "s/^run.*/run $NSTEP/g" \
      step7_production.inp > inp
  sed '/# Running equilibration steps/,+8d' README | \
   sed "s/namd2/\${NAMD} +p1 +devices 0/g" | \
    sed "s/\${prod_prefix}.inp/inp/g" | \
     sed "s/^set cnt    = .*/set cnt = $PREFIX/g" | \
      sed "s/^set cntmax = .*/set cntmax = $PREFIX/g" | \
       grep -v "^#"  >> mdnamd.pbs
  sed -i "s/#PBS -l walltime=.*/#PBS -l walltime=256:00:00/g" mdnamd.pbs
else
  sed -e "s/namd2/\${NAMD} +p4 +devices 0/g" README | grep -v "^#"  >> mdnamd.pbs
fi
