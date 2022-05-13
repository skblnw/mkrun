#!/bin/bash
################################################
# Aug 2015
# CITYUHK-OSU-ZJU
# Kevin
# This script was created for preparing NAMD MD jobs ***IN A LAZY WAY***
#
# Updates:
# May 2022
#  Reduced the script to create generic MD parameters
# May 2017
#  Revised NAMD CUDA-smp-ib commands (highly recommended in the release note)
#  Deleted NVT equilibirum before releasing constraints
#  Added releasing constraints to lipid heads
#  Deleted redundant minimization steps
# 17 Aug 2016
#  Created mknamd_prepare
#  Added function make_pbs
#  Added clusters as options
# 13 Sep 2015
#  Added option bmk
# 4 Aug 2015
#  Added variable previous
# 3 Aug 2015
#  Added variable COMMAND
#  Added variable jobname
#  Added variable nnode and nheader
################################################
# NAMD Commands (GPU):
#   COMBO (remote-shell=ssh): \$CHARMRUN +p\${NP_2} ++ppn $ncpu_2 ++scalable-start ++verbose ++nodelist \$nodefile \$NAMD +idlepoll +setcpuaffinity +pemap 1-$ncpu_2 +commap 0
#   COMBO (remote-shell=ib):  \$CHARMRUN +p\${NP_2} ++ppn $ncpu_2 ++scalable-start ++verbose ++mpiexec ++remote-shell ./mpiexec_wrapper \$NAMD +idlepoll +setcpuaffinity +pemap 1-$ncpu_2 +commap 0
#   Titan: aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0
#   College: mpirun -n \$NP \$NAMD
################################################
# NAMD Commands (CPU):
#   COMBO (remote-shell=ssh): \$CHARMRUN +p\${NP_1} ++ppn $ncpu_1 ++scalable-start ++verbose ++nodelist \$nodefile \$NAMD
#   COMBO (remote-shell=ib):  \$CHARMRUN +p\${NP_1} ++scalable-start ++verbose ++mpiexec ++remote-shell ./mpiexec_wrapper \$NAMD
################################################
# Gromacs Commands (GPU):
#   COMBO: $MPIRUN -mca btl self,openib -np \${NN} -npernode 1 $GMXBIN mdrun -v -ntomp 7 -pin on -s W50k.tpr -deffnm output/gpu-1node-NNomp
################################################


[ $# -ne 2 ] && { echo "mknamd> Usage: $0 [start] [end]"; exit 1; }

START=$1
END=$2

mkdir -p run
for ii in $(seq $START $END); do
    mkdir run/win$ii run/win$ii/pdb2namd run/win$ii/pdb2namd/vmd_solvate
    sed -e 's/CENTER/'$ii'/g' umbrella.col > run/win$ii/umbrella.col
    cd run/win$ii
    ln -s ~/toppar_c36_jul20/ toppar
    ln -s ~/toppar_water_ions_namd.str .
    ln -s ../../template-namd .
    ln -s ../../windows/win$ii.restart.coor smd.restart.coor
    ln -s ../../windows/win$ii.restart.xsc smd.restart.xsc
    cd pdb2namd/vmd_solvate
    ln -s ../../../../ionized.psf .
    ln -s ../../../../ionized.pdb .
    cd ../../../..
done

# /-------------------------------/
# /     Some Global Variables     /
# /-------------------------------/

NORMRUN="namd3 +p2 +devices 0"
CUDARUN="namd3 +p1 +devices 0"
VMD="/opt/vmd/1.9.3/vmd"
RUN_DIR=run

# /-------------------/
# /     Switches      /
# /-------------------/

# MD Steps Choosing
# true to turn on
mini=true
heat=true
cons=true
md=true
md_posres=true
MDPOSRES="segname PROB and resid 23"


# /--------------------------/
# /     Useful Functions     /
# /--------------------------/

# Check existence of a file
check_exist () {
    if [ ! -s "$1" ]; then
        echo "> $1 DOES NOT EXIST! Be careful!"
        exit 1
    fi
}

add_pbs () {
    if [ "$#" -ne 2 ]; then
        echo "> Wrong argument number"
        exit 1
    fi

    echo "echo -e \"mknamd> Working hard with $1...\nmknamd> Finished at:\"" >> $2
    if $acceleration; then
        echo -e "$CUDARUN $1.namd > log/$1.log\nwait\ndate" >> $2
    else
        echo -e "$NORMRUN $1.namd > log/$1.log\nwait\ndate" >> $2
    fi
}

# /---------------------------------/
# /     Prepare restraint files     /
# /---------------------------------/


for windex in $(seq $START $END); do 
    cd run/win${windex}

[ -d "$RUN_DIR" ] && { echo "mknamd> Directory run exists! Backing up to run.BAK"; rm -rf $RUN_DIR.BAK; mv $RUN_DIR $RUN_DIR.BAK; }
[ -d "restraints" ] && { echo "mknamd> Directory restraints exists! Backing up to restraints.BAK"; rm -rf restraints.BAK; mv restraints restraints.BAK; }
mkdir $RUN_DIR $RUN_DIR/output $RUN_DIR/log restraints

cat > tcl <<'EOF'
mol new pdb2namd/vmd_solvate/ionized.psf waitfor all
mol addfile smd.restart.coor waitfor all
set all [atomselect top "all"]

$all set beta 0
set sel [atomselect top {not segname "WT.*" TIP3 ION IONS}]
$sel set beta 1
$all writepdb restraints/cons_solute.pdb

$all set beta 0
set sel [atomselect top {noh and not segname "WT.*" TIP3 ION IONS}]
$sel set beta 1
$all writepdb restraints/cons_heavy.pdb

$all set beta 0
set sel [atomselect top {noh and backbone and not segname "WT.*" TIP3 ION IONS}]
$sel set beta 1
$all writepdb restraints/cons_backbone.pdb

$all set beta 0
set sel [atomselect top {name CA and not segname "WT.*" TIP3 ION IONS}]
$sel set beta 1
$all writepdb restraints/cons_CA.pdb
quit
EOF
$VMD -dispdev text -e tcl > LOG_vmd
rm tcl 


# /-------------------/
# /     Main body     /
# /-------------------/

touch $RUN_DIR/run.sh

# /-----------------------/
# /     Minimization      /
# /-----------------------/
if $mini; then
    acceleration=false
    check_exist template-namd
    prefix=mini
    frequency=5000

    # Change these # of steps. DO NOT continue unless you see a plateau. (Kevin, 2017)
    nstep=(5000 5000 10000)
    ii=1
    for conspdb in cons_solute cons_heavy cons_backbone; do
        jj=$((ii-1))
        if [ $ii -eq 1 ]; then
            inputname="..\/smd.restart"
            outputname="output\/${prefix}-${ii}"
        else
            inputname="output\/${prefix}-${jj}"
            outputname="output\/${prefix}-${ii}"
        fi

        sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
            -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.vel/g' \
            -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.coor/g' \
            -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.xsc/g' \
            -e 's/^set FIXPDB.*$/set FIXPDB 0/g' \
            -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/'${conspdb}'/g' \
            -e 's/^set CONSSCALE.*$/set CONSSCALE 50.0/g' \
            -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
            -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
            -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
            -e 's/^set MS.*$/set MS '${nstep[$jj]}'/g' \
            -e 's/^set SD.*$/set SD 0/g' \
            template-namd > $RUN_DIR/${prefix}-${ii}.namd

        if [ $ii -eq 1 ]; then
            sed -i 's/BinVelocities.*//g' $RUN_DIR/${prefix}-${ii}.namd
        fi

        add_pbs ${prefix}-${ii} $RUN_DIR/run.sh
        ii=$((ii+1))
    done
fi

# /-----------------------/
# /       Annealing       /
# /-----------------------/
if $heat; then
    acceleration=false
    check_exist template-namd
    prefix=heat
    frequency=5000

    previous="output\/mini-3"
    outputname="output\/${prefix}"

    sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${previous}'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
        -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/cons_backbone/g' \
        -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
        -e 's/^set ITEMP.*$/set ITEMP 0/g' \
        -e 's/^set FTEMP.*$/set FTEMP 303/g' \
        -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
        -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
        -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
        -e 's/^set TS.*$/set TS 50000/g' \
        template-namd > $RUN_DIR/${prefix}.namd
    add_pbs ${prefix} $RUN_DIR/run.sh
fi

# /-----------------------------/
# /     Constraining Runs       /
# /-----------------------------/
if $cons; then
    acceleration=true
    check_exist template-namd
    prefix=cons
    frequency=10000

    previous="heat"
    ii=1
    for force in 2 1 .1
    do
        if [ $ii -eq 1 ]; then
            inputname="output\/${previous}"
            outputname="output\/${prefix}-${ii}"
        else
            jj=$((ii-1))
            inputname="output\/${prefix}-${jj}"
            outputname="output\/${prefix}-${ii}"
        fi

        sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
            -e 's/^set CONSSCALE.*$/set CONSSCALE '${force}'/g' \
            -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/cons_CA/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
            -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
            -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
            -e 's/^set TS.*$/set TS 200000/g' \
            -e 's/^set CUDASOA.*$/set CUDASOA 1/g' \
            template-namd > $RUN_DIR/${prefix}-${ii}.namd

        add_pbs ${prefix}-${ii} $RUN_DIR/run.sh
        ii=$((ii+1))
    done
fi

# /-----------------------------/
# /      Umbrella Sampling      /
# /-----------------------------/
if $md; then
    acceleration=false
    check_exist template-namd
    prefix=md
    frequency=10000

    inputname="output\/cons-3"
    outputname="output\/${prefix}"

    sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
        -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.restart.vel/g' \
        -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.restart.coor/g' \
        -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.restart.xsc/g' \
        -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
        -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
        -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
        -e 's/^COMmotion.*$/COMmotion no/g' \
        -e 's/^set ITEMP.*$/set ITEMP 303/g' \
        -e 's/^set FTEMP.*$/set FTEMP 303/g' \
        -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
        -e 's/^set TS.*$/set TS 1000000/g' \
        -e 's/^set CUDASOA.*$/set CUDASOA 0/g' \
        -e 's/^set CONSSCALE.*$/set CONSSCALE 10/g' \
        -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/cons_posres/g' \
        -e 's/^set COLFILE.*$/set COLFILE ..\/umbrella.col/g' \
        template-namd > $RUN_DIR/${prefix}.namd

    cat > tcl <<EOF
mol new pdb2namd/vmd_solvate/ionized.psf waitfor all
mol addfile smd.restart.coor waitfor all
set all [atomselect top "all"]
\$all set beta 0
set sel [atomselect top "$MDPOSRES and name CA"]
\$sel set beta 1
\$all writepdb restraints/cons_posres.pdb
quit
EOF
    $VMD -dispdev text -e tcl >> LOG_vmd; rm tcl

    add_pbs ${prefix} $RUN_DIR/run.sh
fi

cd ../..
done
