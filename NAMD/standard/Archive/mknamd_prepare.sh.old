#!/bin/bash
################################################
# Aug 2015
# CBMSG
# Kevin
# This script was created for preparing NAMD MD jobs ***IN A LAZY WAY***
#
# Updates:
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
# /     RUN Command Switches      /
# /-------------------------------/

NORMRUN="namd3 +p4 +devices 0"
CUDARUN="namd3 +p1 +devices 0"
VMD="/opt/vmd/1.9.3/vmd"

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
MDPOSRES="segname PROA and resid 90"
md_continue=false
bmk=false
smd=false
us1=false
us2=false




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

# Create a pbs acoording to chosen cluster
# Usage: make_pbs $jobname ($nnode_total)
make_pbs () {
    check_exist template-job-pbs

    if [ "$#" -ne 1 ]; then
        echo "> Parameter \$jobname is empty!"
        exit 1
    fi

    sed -e 's/^#PBS -N .*/#PBS -N '$1'/g' template-job-pbs > job-$1.pbs

    if $titan; then
        if [ -z "$2" ]; then
            sed -i 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode'/g' job-$1.pbs
        else
            sed -i 's/^#PBS -l nodes=.*/#PBS -l nodes='$2'/g' job-$1.pbs
        fi
    else
        sed -i 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' job-$1.pbs
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

[ -d "run" ] && { echo "mknamd> Directory run exists! Backing up to run.BAK"; rm -rf run.BAK; mv run run.BAK; }
[ -d "restraints" ] && { echo "mknamd> Directory restraints exists! Backing up to restraints.BAK"; rm -rf restraints.BAK; mv restraints restraints.BAK; }
mkdir run run/output run/log restraints
touch run/run.sh

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

# /--------------------------------------/
# /     Only if you run on a cluster     /
# /--------------------------------------/

run_on_cluster=false
nnode=2
ncpu=10
ncpu_1=$(($ncpu-1))
ncpu_2=$(($ncpu-2))

# Cluster Choosing
# true if you choose that cluster
combo=false
# Combo has a cpu version of pbs
combo_cpu=false
titan=false
# Titan has a multiple version of pbs
titan_multiple=false
college=false
chenglab=true
if $combo; then
    if [[ $nnode -eq 1 ]]; then
        COMMAND="\$NAMD +p$ncpu_1"
    elif $combo_cpu; then
        COMMAND="\$CHARMRUN +p$ncpu_1 ++scalable-start ++verbose ++mpiexec ++remote-shell ./mpiexec_wrapper \$NAMD"
    else
        COMMAND="\$CHARMRUN +p\${NP_2} ++ppn $ncpu_2 ++scalable-start ++verbose ++mpiexec ++remote-shell ./mpiexec_wrapper \$NAMD +idlepoll +setcpuaffinity +pemap 1-$ncpu_2 +commap 0"
    fi
elif $titan; then
    COMMAND="aprun -n $nnode -r 1 -N 1 -d 15 namd2 +ppn 15 +pemap 1-15 +commap 0"
elif $college; then
    COMMAND="mpirun -n \$NP \$NAMD"
elif $chenglab; then
    COMMAND="\$NAMD +p$ncpu"
else
    echo "Specify on which machine you are gonna run!"
    exit 0
fi

COMMAND="namd3 +p1 +devices 0"

# /-----------------------/
# /     Minimization      /
# /-----------------------/
if $mini; then
    acceleration=false
    check_exist template-namd
    prefix=mini
    frequency=5000

    if $run_on_cluster; then
        make_pbs $prefix
    fi

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
            template-namd > run/${prefix}-${ii}.namd

        add_pbs ${prefix}-${ii} run/run.sh
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

    if $run_on_cluster; then
        make_pbs $prefix
    else 
        touch run/run.sh
    fi

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
        template-namd > run/${prefix}.namd
    add_pbs ${prefix} run/run.sh
fi

# /-----------------------------/
# /     Constraining Runs       /
# /-----------------------------/
if $cons; then
    acceleration=false
    check_exist template-namd
    prefix=cons
    frequency=10000

    if $run_on_cluster; then
        make_pbs $jobname
    fi

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
                template-namd > run/${prefix}-${ii}.namd

            add_pbs ${prefix}-${ii} run/run.sh
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
                template-namd > run/${prefix}-${ii}.namd

            add_pbs ${prefix}-${ii} run/run.sh
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

    if $run_on_cluster; then
        make_pbs $jobname
    fi

    inputname="output\/cons-6"
    outputname="output\/${prefix}"

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
        sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
            -e 's/^set CONSSCALE.*$/set CONSSCALE 10/g' \
            -e 's/^set CONSPDB.*$/set CONSPDB ..\/restraints\/cons_posres/g' \
            -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.restart.vel/g' \
            -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.restart.coor/g' \
            -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.restart.xsc/g' \
            -e 's/^COMmotion.*$/COMmotion no/g' \
            -e 's/^set ITEMP.*$/set ITEMP 303/g' \
            -e 's/^set FTEMP.*$/set FTEMP 303/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
            -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
            -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
            -e 's/^set TS.*$/set TS 50000000/g' \
            -e 's/^set CUDASOA.*$/set CUDASOA 1/g' \
            template-namd > run/${prefix}.namd
        add_pbs ${prefix} run/run.sh
    else
        sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
            -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.restart.vel/g' \
            -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.restart.coor/g' \
            -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.restart.xsc/g' \
            -e 's/^COMmotion.*$/COMmotion no/g' \
            -e 's/^set ITEMP.*$/set ITEMP 303/g' \
            -e 's/^set FTEMP.*$/set FTEMP 303/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
            -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
            -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
            -e 's/^set TS.*$/set TS 50000000/g' \
            -e 's/^set CUDASOA.*$/set CUDASOA 1/g' \
            template-namd > run/${prefix}.namd
        add_pbs ${prefix} run/run.sh
    fi
fi

if $md_continue; then
    acceleration=true
    check_exist template-namd
    prefix=md
    frequency=50000
    for ii in {2..3}
    do

        jj=$((ii-1))
        inputname="output\/${prefix}-${jj}"
        outputname="output\/${prefix}-${ii}"

        sed \
            -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME '${outputname}'/g' \
            -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.restart.vel/g' \
            -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.restart.coor/g' \
            -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.restart.xsc/g' \
            -e 's/^COMmotion.*$/COMmotion yes/g' \
            -e 's/^set ITEMP.*$/set ITEMP 310/g' \
            -e 's/^set FTEMP.*$/set FTEMP 310/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq '${frequency}'/g' \
            -e 's/^dcdfreq.*$/dcdfreq '${frequency}'/g' \
            -e 's/^xstfreq.*$/xstfreq '${frequency}'/g' \
            -e 's/^set TS.*$/set TS 25000000/g' \
            template-namd > ${prefix}-$ii.namd
        add_pbs ${prefix}-${ii} run/run.sh

    done
fi

if $bmk; then
    check_exist template-namd
    ncpu=16
    for nnode in {1,2,8,16,32,64}
    do
        jobname=bmk$nnode
        make_pbs $jobname

        previous=NPT-1
        sed \
            -e 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/bmk-'$nnode'/g' \
            -e 's/^COMmotion.*$/COMmotion yes/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq 1000/g' \
            -e 's/^dcdfreq.*$/dcdfreq 1000/g' \
            -e 's/^xstfreq.*$/xstfreq 1000/g' \
            -e 's/^set TS.*$/set TS 5000/g' \
            template-namd > bmk-$nnode.namd
        echo -e "$COMMAND bmk-$nnode.namd > bmk-$nnode.log\nwait\ndate" >> job-$nnode.pbs
    done
fi

if $smd; then
    check_exist template-smd-namd

    # set total SMD distance (A)
    dist=20
    # set SMD speed (A/ns)
    for spd in
    do
        filename=smd$spd
    
        previous=cons-1
        for ii in {1..5}
        do
            jobname=smd$spd-$ii
            make_pbs $jobname
            
            jj=$(expr $ii - 1)
            sed -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/'$filename'-'$ii'/g' \
                -e 's/^set FIXPDB.*$/set FIXPDB restraints\/fix_memb/g' \
                -e 's/^COMmotion.*$/COMmotion no/g' \
                -e 's/^set PSWITCH.*$/set PSWITCH 0/g' \
                -e 's/^restartfreq.*$/restartfreq 1000/g' \
                -e 's/^dcdfreq.*$/dcdfreq 1000/g' \
                -e 's/^xstfreq.*$/xstfreq 1000/g' \
                -e 's/^SMDOutputFreq.*$/SMDOutputFreq 1000/g' \
                -e 's/set SEP.*$/set SEP '$dist'/g' \
                -e 's/set SPEED_GLOB.*$/set SPEED_GLOB '$spd'/g' \
                template-smd-namd > $filename-$ii.namd
                
            if [ $ii -eq 1 ] ; then
                sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' $filename-$ii.namd
            else
                sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/'$filename'-'$jj'/g' $filename-$ii.namd
            fi
        
            echo -e "$COMMAND $filename-$ii.namd > log/$filename-$ii.log\nwait\ndate" >> job-$jobname.pbs
        done
    done
fi

if $us1; then
    check_exist template-us-namd

    if $titan_multiple; then
        read -p "How many node in total will be used? `echo $'\n> '`" nnode_total
        if ! [[ $nnode_total =~ ^-?[0-9]+$ ]]; then
            echo "Must be an integer!"
            exit 1
        fi

        jobname=us1
        make_pbs $jobname $nnode_total
    fi

    for window in 
    do
        OUTPUT_DIR=us-z$window
        mkdir -p $OUTPUT_DIR
        filename=us-z$window

        center_aa=`echo "scale=1; $window / 10" | bc`
        sed -e 's/CENTER/'$center_aa'/g' \
            template-colvars > $OUTPUT_DIR/$filename.col

        sed -e 's/^set INPUTNAME.*$/set INPUTNAME 0/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME ..\/output\/'$filename'-1/g' \
            -e 's/BinCoordinates/BinCoordinates ..\/output\/'$filename'-0.restart.coor/g' \
            -e 's/^colvarsConfig.*$/colvarsConfig us-z'$window'.col/g' \
            -e 's/^COMmotion.*$/COMmotion no/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq 1000/g' \
            -e 's/^dcdfreq.*$/dcdfreq 1000/g' \
            -e 's/^xstfreq.*$/xstfreq 1000/g' \
            -e 's/^set TS.*$/set TS 500000/g' \
            template-us-namd > $OUTPUT_DIR/$filename-1.namd
        if $titan_multiple; then
            echo -e "$COMMAND $OUTPUT_DIR/$filename-1.namd > log/$filename-1.log &\n" >> job-$jobname.pbs
        else
            jobname=z$window-1
            make_pbs $jobname
            echo -e "$COMMAND $OUTPUT_DIR/$filename-1.namd > log/$filename-1.log\nwait\ndate" >> job-$jobname.pbs
        fi
    done

    if $titan_multiple; then
        echo -e "\nwait\ndate" >> job-$jobname.pbs
    fi
fi

if $us2; then
    check_exist template-us-namd

    if $titan_multiple; then
        read -p "How many node in total will be used? `echo $'\n> '`" nnode_total
        if ! [[ $nnode_total =~ ^-?[0-9]+$ ]]; then
            echo "Must be an integer!"
            exit 1
        fi
    fi

    for ii in {2..12}
    do
        jj=$(expr $ii - 1)
        if $titan_multiple; then
            jobname=us$ii
            make_pbs $jobname $nnode_total
        fi
    
        for window in 
        do
            OUTPUT_DIR=us-z$window
            mkdir -p $OUTPUT_DIR
            filename=us-z$window
    
            sed -e 's/^set INPUTNAME.*$/set INPUTNAME ..\/output\/'$filename'-'$jj'/g' \
                -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME ..\/output\/'$filename'-'$ii'/g' \
                -e 's/^colvarsConfig.*$/colvarsConfig us-z'$window'.col/g' \
                -e 's/^COMmotion.*$/COMmotion yes/g' \
                -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
                -e 's/^restartfreq.*$/restartfreq 1000/g' \
                -e 's/^dcdfreq.*$/dcdfreq 1000/g' \
                -e 's/^xstfreq.*$/xstfreq 1000/g' \
                -e 's/^set TS.*$/set TS 5000000/g' \
                template-us-namd > $OUTPUT_DIR/$filename-$ii.namd
            if $titan_multiple; then
                echo -e "$COMMAND $OUTPUT_DIR/$filename-$ii.namd > log/$filename-$ii.log &" >> job-$jobname.pbs
            else
                jobname=z$window-$ii
                make_pbs $jobname
                echo -e "$COMMAND $OUTPUT_DIR/$filename-$ii.namd > log/$filename-$ii.log\nwait\ndate" >> job-$jobname.pbs
            fi
        done
    
        if $titan_multiple; then
            echo -e "\nwait\ndate" >> job-$jobname.pbs
        fi
    done
fi
