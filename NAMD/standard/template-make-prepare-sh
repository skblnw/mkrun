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

mkdir -p output
mkdir -p log
mkdir -p restraints

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


membrane_exist=false
# MD Steps Choosing
# true to turn on
mini=false
heat=false
pre=false
cons=false
npt1=false
npt2=false
bmk=false
smd=false
us1=false
us2=false

# Check existence of a file
check_exist () {
    if [ ! -s "$1" ]; then
        echo "$1 DOES NOT EXIST! Be careful!"
        exit 1
    fi
}

# Create a pbs acoording to chosen cluster
# Usage: make_pbs $jobname ($nnode_total)
make_pbs () {
    check_exist template-job-pbs

    if [ -z "$1" ]; then
        echo "-Parameter \$jobname is empty!"
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


# Minimization
if $mini; then
    check_exist template-namd
    jobname=mini
    make_pbs $jobname
    outputname=$jobname

    fixpdb=(fix_solute fix_heavy 0)
    # Change these # of steps. DO NOT continue unless you see a plateau. (Kevin, 2017)
    nstep=(5000 20000 50000)
    for ii in {0..2}
    do
        if [ $ii -eq 0 ]; then
            inputname=0
        else
            jj=$[ii-1]
            inputname="output\/${outputname}-${jj}"
        fi
        sed -e 's/^set INPUTNAME.*$/set INPUTNAME '${inputname}'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/'${outputname}'-'${ii}'/g' \
            -e 's/^set FIXPDB.*$/set FIXPDB '${fixpdb[$ii]}'/g' \
            -e 's/^set CONSPDB.*$/set CONSPDB restraints\/cons_bb_and_P/g' \
            -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
            -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.vel/g' \
            -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.coor/g' \
            -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.xsc/g' \
            -e 's/^restartfreq.*$/restartfreq 5000/g' \
            -e 's/^dcdfreq.*$/dcdfreq 5000/g' \
            -e 's/^xstfreq.*$/xstfreq 5000/g' \
            -e 's/^set MS.*$/set MS '${nstep[$ii]}'/g' \
            -e 's/^set SD.*$/set SD 0/g' \
            template-namd > ${outputname}-${ii}.namd
        echo -e "$COMMAND ${outputname}-${ii}.namd > log/${outputname}-${ii}.log\nwait\ndate" >> job-$jobname.pbs
    done
fi

if $heat; then
    check_exist template-namd
    jobname=heat
    make_pbs $jobname
    outputname=$jobname

    previous=mini-2
    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/'${outputname}'/g' \
        -e 's/^set CONSPDB.*$/set CONSPDB restraints\/cons_bb_and_P/g' \
        -e 's/^set CONSSCALE.*$/set CONSSCALE 10.0/g' \
        -e 's/^set ITEMP.*$/set ITEMP 0/g' \
        -e 's/^set FTEMP.*$/set FTEMP 310/g' \
        -e 's/^set TS.*$/set TS 100000/g' \
        template-namd > ${outputname}.namd
    echo -e "$COMMAND ${outputname}.namd > log/${outputname}.log\nwait\ndate" >> job-$jobname.pbs
fi

if $pre; then
    check_exist template-namd
    jobname=pre
    make_pbs $jobname

    outputname=$jobname-npt
    previous=heat
    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/'${outputname}'/g' \
        -e 's/^set CONSPDB.*$/set CONSPDB restraints\/cons_bb_and_P/g' \
        -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
        -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
        -e 's/^margin.*$/margin 5.0/g' \
        -e 's/^set TS.*$/set TS 500000/g' \
        template-namd > ${outputname}.namd
    echo -e "$COMMAND ${outputname}.namd > log/${outputname}.log\nwait\ndate" >> job-$jobname.pbs

#    for ii in {1..5}
#    do
#        jj=$(expr $ii - 1)
#        sed -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/pre-npt-cont'$ii'/g' \
#            -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
#            -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
#            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
#            -e 's/^set TS.*$/set TS 50000/g' \
#            template-namd > pre-npt-cont$ii.namd
#        if [ $jj -eq 0 ] ; then
#            sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/pre-npt/g' pre-npt-cont$ii.namd
#        else
#            sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/pre-npt-cont'$jj'/g' pre-npt-cont$ii.namd
#        fi
#        echo -e "$COMMAND pre-npt-cont$ii.namd > log/pre-npt-cont$ii.log\nwait\ndate" >> job-$jobname.pbs
#    done

#    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/pre-npt-cont5/g' \
#        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/pre-nvt/g' \
#        -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
#        -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
#        -e 's/^set PSWITCH.*$/set PSWITCH 0/g' \
#        -e 's/^set TS.*$/set TS 50000/g' \
#        template-pre > pre-nvt.namd
#    echo -e "$COMMAND pre-nvt.namd > log/pre-nvt.log\nwait\ndate" >> job-$jobname.pbs
fi

if $cons; then
    check_exist template-namd
    jobname=cons
    make_pbs $jobname

    MM=-1
    previous=pre-npt

    if $membrane_exist; then
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
    fi

    outputname=$jobname
    for NN in {50,20,10,5,2,1}
    do
        SS=`echo "$NN 0.1" | awk '{printf "%.1f", $1*$2}'`
        sed   -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/'${outputname}'-'$NN'/g' \
              -e 's/^set CONSSCALE.*$/set CONSSCALE '$SS'/g' \
              -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
              -e 's/^set TS.*$/set TS 50000/g' \
              template-namd > cons-$NN.namd

        if [ $MM -eq -1 ] ; then
            sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' ${outputname}-$NN.namd
        else
            sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/'${outputname}'-'$MM'/g' ${outputname}-$NN.namd
        fi

        if [ $NN -gt 5 ] ; then
            sed -i 's/^set CONSPDB.*$/set CONSPDB restraints\/cons_bb/g' ${outputname}-$NN.namd
        else
            sed -i 's/^set CONSPDB.*$/set CONSPDB restraints\/cons_CA/g' ${outputname}-$NN.namd
        fi

        echo -e "$COMMAND ${outputname}-$NN.namd > log/${outputname}-$NN.log\nwait\ndate" >> job-$jobname.pbs
        MM=$NN
    done
fi

if $npt1; then
    check_exist template-namd
    jobname=npt1
    make_pbs $jobname
    outputname=NPT

    previous=cons-1
    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/'${outputname}'-1/g' \
        -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.restart.vel.old/g' \
        -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.restart.coor.old/g' \
        -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.restart.xsc.old/g' \
        -e 's/^COMmotion.*$/COMmotion no/g' \
        -e 's/^set ITEMP.*$/set ITEMP 310/g' \
        -e 's/^set FTEMP.*$/set FTEMP 310/g' \
        -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
        -e 's/^restartfreq.*$/restartfreq 1000/g' \
        -e 's/^dcdfreq.*$/dcdfreq 1000/g' \
        -e 's/^xstfreq.*$/xstfreq 1000/g' \
        -e 's/^set TS.*$/set TS 500000/g' \
        template-namd > ${outputname}-1.namd
    echo -e "$COMMAND ${outputname}-1.namd > log/${outputname}-1.log\nwait\ndate" >> job-$jobname.pbs
fi

if $npt2; then
    check_exist template-namd
    for ii in {2..3}
    do
        jobname=npt$ii
        make_pbs $jobname
        outputname=NPT

        jj=$(expr $ii - 1)
        sed \
            -e 's/^set INPUTNAME.*$/set INPUTNAME output\/'${outputname}'-'$jj'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/'${outputname}'-'$ii'/g' \
            -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.restart.vel.old/g' \
            -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.restart.coor.old/g' \
            -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.restart.xsc.old/g' \
            -e 's/^COMmotion.*$/COMmotion yes/g' \
            -e 's/^set ITEMP.*$/set ITEMP 310/g' \
            -e 's/^set FTEMP.*$/set FTEMP 310/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq 1000/g' \
            -e 's/^dcdfreq.*$/dcdfreq 1000/g' \
            -e 's/^xstfreq.*$/xstfreq 1000/g' \
            -e 's/^set TS.*$/set TS 25000000/g' \
            template-namd > ${outputname}-$ii.namd
        echo -e "$COMMAND ${outputname}-$ii.namd > log/${outputname}-$ii.log\nwait\ndate" >> job-$jobname.pbs
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
