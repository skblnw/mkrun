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




# /-------------------------------/
# /     Some Global Variables     /
# /-------------------------------/

NORMRUN="namd3 +p4 +devices 0"
CUDARUN="namd3 +p1 +devices 0"
VMD="/opt/vmd/1.9.3/vmd"
RUN_DIR=run

# /-------------------/
# /     Switches      /
# /-------------------/

membrane_exist=false
# MD Steps Choosing
# true to turn on
mini=true
heat=true
cons=true
md=true
md_posres=true
MDPOSRES="segname PROA PROB"
smd=false
SMDGROUP="segname PROC PROD and name CA"


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

[ -d "$RUN_DIR" ] && { echo "mknamd> Directory run exists! Backing up to run.BAK"; rm -rf $RUN_DIR.BAK; mv $RUN_DIR $RUN_DIR.BAK; }
[ -d "restraints" ] && { echo "mknamd> Directory restraints exists! Backing up to restraints.BAK"; rm -rf restraints.BAK; mv restraints restraints.BAK; }
mkdir $RUN_DIR $RUN_DIR/output $RUN_DIR/log restraints

cat > tcl <<'EOF'
mol new pdb2namd/vmd_solvate/ionized.pdb type pdb waitfor all
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

COMMAND="namd3 +p1 +devices 0"
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
    for conspdb in cons_solute cons_heavy cons_backbone
    do
        jj=$((ii-1))
        if [ $ii -eq 1 ]; then
            inputname=0
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

    if $membrane_exist; then

        cat > tcl <<EOF
mol new pdb2namd/vmd_solvate/ionized.pdb type pdb waitfor all
set all [atomselect top "all"]

\$all set beta 0
set sel [atomselect top "segname MEMB and name P"]
\$sel get resid
\$sel get resname
\$sel get name
\$sel set beta 1
\$all writepdb restraints/cons_posres.pdb
quit
EOF

        outputname=$jobname-memb
        for NN in {5,2,1}
        do
            # Get membrane_lipid_restraint.namd.col from CHARMM-GUI
            colfile=${outputname}-$NN.col
            sed -e 's/Constant \$fc/Constant '$NN'/g' membrane_lipid_restraint.namd.col > restraints/$colfile

            sed   -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/'${outputname}'-'$NN'/g' \
                  -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
                  -e 's/^set CONSPDB.*$/set CONSPDB restraints\/cons_bb/g' \
                  -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
                  -e 's/^set COLFILE.*$/set COLFILE restraints\/'$colfile'/g' \
                  template-namd > ${outputname}-$NN.namd

            if [ $NN -eq 5 ]; then
                sed -i 's/^set TS.*$/set TS 100000/g' ${outputname}-$NN.namd
            else
                sed -i 's/^set TS.*$/set TS 25000/g' ${outputname}-$NN.namd
            fi

            if [ $MM -eq -1 ] ; then
                sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' ${outputname}-$NN.namd
            else
                sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/'${outputname}'-'$MM'/g' ${outputname}-$NN.namd
            fi

            echo -e "$COMMAND ${outputname}-$NN.namd > log/${outputname}-$NN.log\nwait\ndate" >> job-$jobname.pbs
            MM=$NN
        done
        MM=-1
        previous=${outputname}-${NN}

    else  

        previous="heat"
        ii=1
        for force in 5 2 1
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
                -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/cons_backbone/g' \
                -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
                -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
                -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
                -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
                -e 's/^set TS.*$/set TS 100000/g' \
                -e 's/^set CUDASOA.*$/set CUDASOA 1/g' \
                template-namd > $RUN_DIR/${prefix}-${ii}.namd

            add_pbs ${prefix}-${ii} $RUN_DIR/run.sh
            ii=$((ii+1))
        done

        for force in 1 .5 .1
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
                -e 's/^set TS.*$/set TS 100000/g' \
                -e 's/^set CUDASOA.*$/set CUDASOA 1/g' \
                template-namd > $RUN_DIR/${prefix}-${ii}.namd

            add_pbs ${prefix}-${ii} $RUN_DIR/run.sh
            ii=$((ii+1))
        done

    fi
fi

# /---------------------------/
# /     Production Runs       /
# /---------------------------/
if $md; then
    acceleration=true
    check_exist template-namd
    prefix=md
    frequency=50000

    inputname="output\/cons-6"
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
        -e 's/^set TS.*$/set TS 50000000/g' \
        -e 's/^set CUDASOA.*$/set CUDASOA 1/g' \
        template-namd > namdpar

    if $md_posres; then
        cat > tcl <<EOF
mol new pdb2namd/vmd_solvate/ionized.pdb type pdb waitfor all
set all [atomselect top "all"]

\$all set beta 0
set sel [atomselect top "$MDPOSRES and name CA"]
\$sel get resid
\$sel get resname
\$sel get name
\$sel set beta 1
\$all writepdb restraints/cons_posres.pdb
quit
EOF
        $VMD -dispdev text -e tcl >> LOG_vmd; rm tcl

        sed -e 's/^set CONSSCALE.*$/set CONSSCALE 10/g' \
            -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/cons_posres/g' \
            namdpar > namdpar2
        mv namdpar2 namdpar
    fi
    mv namdpar $RUN_DIR/${prefix}.namd
    add_pbs ${prefix} $RUN_DIR/run.sh
fi

# /---------------------------/
# /        Steered MD         /
# /---------------------------/
if $smd; then
    acceleration=false
    check_exist template-namd
    check_exist smd.col
    prefix=smd
    frequency=10000

    inputname="output\/cons-6"
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
        -e 's/^set TS.*$/set TS 300000/g' \
        -e 's/^set CONSSCALE.*$/set CONSSCALE 10/g' \
        -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/cons_posres/g' \
        -e 's/^set SMDFILE.*$/set SMDFILE ..\/restraints\/ionized.smd/g' \
        -e 's/^set COLFILE.*$/set COLFILE ..\/smd.col/g' \
        template-namd > $RUN_DIR/${prefix}.namd

    cat > tcl <<EOF
mol new pdb2namd/vmd_solvate/ionized.pdb type pdb waitfor all
set all [atomselect top "all"]

\$all set beta 0
set sel [atomselect top "$MDPOSRES and name CA"]
\$sel set beta 1
\$all writepdb restraints/cons_posres.pdb

\$all set occupancy 0
set sel [atomselect top "$SMDGROUP"]
\$sel set occupancy 1
\$all writepdb restraints/ionized.smd
quit
EOF
    $VMD -dispdev text -e tcl >> LOG_vmd; rm tcl

    add_pbs ${prefix} $RUN_DIR/run.sh
fi
