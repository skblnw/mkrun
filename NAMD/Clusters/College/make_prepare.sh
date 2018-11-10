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

nnode=5
ncpu=16
nheader=30
COMMAND="\$MPIRUN -n \$NP \$NAMD"

mini=true
heat=false
pre=false
cons=false
npt1=false
npt2=false
bmk=false

if $mini; then
    jobname=mini
    sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
        -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
        -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
        template-job-pbs > job-$jobname.pbs
    sed -i $nheader',$d' job-$jobname.pbs

    fixpdb=(fix_solute 0)
    nstep=(1000 20000)
    for ii in {0..1}
    do
        if [ $ii -eq 0 ]; then
            inputname=0
        else
            jj=$[ii-1]
            inputname="output\/mini-$jj"
        fi
        sed -e 's/^set INPUTNAME.*$/set INPUTNAME '$inputname'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/mini-'$ii'/g' \
            -e 's/^set FIXPDB.*$/set FIXPDB '${fixpdb[$ii]}'/g' \
            -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
            -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
            -e 's/BinVelocities.*$/BinVelocities $INPUTNAME.vel/g' \
            -e 's/BinCoordinates.*$/BinCoordinates $INPUTNAME.coor/g' \
            -e 's/ExtendedSystem.*$/ExtendedSystem $INPUTNAME.xsc/g' \
            -e 's/^restartfreq.*$/restartfreq 1000/g' \
            -e 's/^dcdfreq.*$/dcdfreq 1000/g' \
            -e 's/^xstfreq.*$/xstfreq 1000/g'
            -e 's/^set MS.*$/set MS '${nstep[$ii]}'/g' \
            -e 's/^set SD.*$/set SD 0/g' \
            template-namd > mini-$ii.namd
        echo -e "$COMMAND mini-$ii.namd > mini-$ii.log\nwait" >> job-$jobname.pbs
    done

    # NN=1
    # fixname=fix_solute
    # sed -e 's/^set INPUTNAME.*$/set INPUTNAME 0/g' \
    #     -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/mini-'$NN'/g' \
    #     -e 's/^set FIXPDB.*$/set FIXPDB '$fixname'/g' \
    #     -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
    #     -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
    #     -e 's/^set MS.*$/set MS 2000/g' \
    #     -e 's/^set SD.*$/set SD 0/g' \
    #     template-namd > mini-$NN.namd
    # echo -e "$COMMAND mini-$NN.namd > mini-$NN.log\nwait" >> job-$jobname.pbs

    # NN=2
    # fixname=fix_all0
    # sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/mini-1/g' \
    #     -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/mini-'$NN'/g' \
    #     -e 's/^set FIXPDB.*$/set FIXPDB 0/g' \
    #     -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
    #     -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
    #     -e 's/^set MS.*$/set MS 2000/g' \
    #     -e 's/^set SD.*$/set SD 0/g' \
    #     template-namd > mini-$NN.namd
    # echo -e "$COMMAND mini-$NN.namd > mini-$NN.log\nwait" >> job-$jobname.pbs
#
#    NN=2
#    MM=$(expr $NN - 1)
#    fixname=fix_backbone
#    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/mini-'$MM'/g' \
#        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/mini-'$NN'/g' \
#        -e 's/^set FIXPDB.*$/set FIXPDB '$fixname'/g' \
#        -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
#        -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
#        -e 's/^set MS.*$/set MS 5000/g' \
#        -e 's/^set SD.*$/set SD 0/g' \
#        template-namd > mini-$NN.namd
#    echo -e "\$NAMD +ppn $ncpu +pemap 1-$ncpu +commap 0 mini-$NN.namd > mini-$NN.log\nwait" >> job-mini.pbs
#
#    NN=3
#    MM=$(expr $NN - 1)
#    fixname=fix_water
#    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/mini-'$MM'/g' \
#        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/mini-'$NN'/g' \
#        -e 's/^set FIXPDB.*$/set FIXPDB '$fixname'/g' \
#        -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
#        -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
#        -e 's/^set MS.*$/set MS 10000/g' \
#        -e 's/^set SD.*$/set SD 0/g' \
#        template-namd > mini-$NN.namd
#    echo -e "\$NAMD +ppn $ncpu +pemap 1-$ncpu +commap 0 mini-$NN.namd > mini-$NN.log\nwait" >> job-mini.pbs
fi

if $heat; then
    jobname=heat
    sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
        -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
        -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
        template-job-pbs > job-$jobname.pbs
    sed -i $nheader',$d' job-$jobname.pbs

    previous=mini-1
    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/heat/g' \
        -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
        -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
        -e 's/^set ITEMP.*$/set ITEMP 0/g' \
        -e 's/^set TS.*$/set TS 50000/g' \
        template-namd > heat.namd
    echo -e "$COMMAND heat.namd > heat.log\nwait" >> job-$jobname.pbs
fi

if $pre; then
    jobname=pre
    sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
        -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
        -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
        template-job-pbs > job-$jobname.pbs
    sed -i $nheader',$d' job-$jobname.pbs

    previous=heat
    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/pre-npt/g' \
        -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
        -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
        -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
        -e 's/^margin.*$/margin 3.0/g' \
        -e 's/^set TS.*$/set TS 500000/g' \
        template-namd > pre-npt.namd
    echo -e "$COMMAND pre-npt.namd > pre-npt.log\nwait" >> job-$jobname.pbs

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
#        echo -e "$COMMAND pre-npt-cont$ii.namd > pre-npt-cont$ii.log\nwait" >> job-pre.pbs
#    done

#    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/pre-npt-cont5/g' \
#        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/pre-nvt/g' \
#        -e 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' \
#        -e 's/^set CONSSCALE.*$/set CONSSCALE 5.0/g' \
#        -e 's/^set PSWITCH.*$/set PSWITCH 0/g' \
#        -e 's/^set TS.*$/set TS 50000/g' \
#        template-pre > pre-nvt.namd
#    echo -e "$COMMAND pre-nvt.namd > pre-nvt.log\nwait" >> job-pre.pbs
fi

MM=-1
if $cons; then
    jobname=cons
    sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
        -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
        -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
        template-job-pbs > job-$jobname.pbs
    sed -i $nheader',$d' job-$jobname.pbs

    previous=pre-npt
    for NN in {50,30,10,5,3,1}
    do
        SS=`echo "$NN 0.1" | awk '{printf "%.1f", $1*$2}'`
        sed   -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/cons-'$NN'/g' \
              -e 's/^set CONSSCALE.*$/set CONSSCALE '$SS'/g' \
              -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
              -e 's/^set TS.*$/set TS 50000/g' \
              template-namd > cons-$NN.namd

        if [ $MM -eq -1 ] ; then
            sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' cons-$NN.namd
        else
            sed -i 's/^set INPUTNAME.*$/set INPUTNAME output\/cons-'$MM'/g' cons-$NN.namd
        fi

        if [ $NN -gt 5 ] ; then
            sed -i 's/^set CONSPDB.*$/set CONSPDB cons_bb_and_P/g' cons-$NN.namd
        else
            sed -i 's/^set CONSPDB.*$/set CONSPDB cons_CA_and_P/g' cons-$NN.namd
        fi

        echo -e "$COMMAND cons-$NN.namd > cons-$NN.log\nwait" >> job-$jobname.pbs
        MM=$NN
    done
fi

if $npt1; then
    jobname=npt-1
    sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
        -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
        -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
        template-job-pbs > job-1.pbs
    sed -i $nheader',$d' job-1.pbs

    previous=cons-1
    sed -e 's/^set INPUTNAME.*$/set INPUTNAME output\/'$previous'/g' \
        -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/NPT-1/g' \
        -e 's/^COMmotion.*$/COMmotion no/g' \
        -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
        -e 's/^set TS.*$/set TS 500000/g' \
        template-namd > NPT-1.namd
    echo -e "$COMMAND NPT-1.namd > NPT-1.log\nwait" >> job-1.pbs
fi

if $npt2; then
    for ii in {2..5}
    do
        jobname=npt-$ii
        sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
            -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
            -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=02:00:00/g' \
            template-job-pbs > job-$ii.pbs
        sed -i $nheader',$d' job-$ii.pbs

        jj=$(expr $ii - 1)
        sed \
            -e 's/^set INPUTNAME.*$/set INPUTNAME output\/NPT-'$jj'/g' \
            -e 's/^set OUTPUTNAME.*$/set OUTPUTNAME output\/NPT-'$ii'/g' \
            -e 's/^COMmotion.*$/COMmotion yes/g' \
            -e 's/^set PSWITCH.*$/set PSWITCH 1/g' \
            -e 's/^restartfreq.*$/restartfreq 1000/g' \
            -e 's/^dcdfreq.*$/dcdfreq 1000/g' \
            -e 's/^xstfreq.*$/xstfreq 1000/g' \
            -e 's/^set TS.*$/set TS 5000000/g' \
            template-namd > NPT-$ii.namd
        echo -e "$COMMAND NPT-$ii.namd > NPT-$ii.log\nwait" >> job-$ii.pbs
    done
fi

if $bmk; then
    ncpu=16
    for nnode in {1,2,8,16,32,64}
    do
        jobname=bmk-$nnode
        sed -e 's/^#PBS -N .*/#PBS -N '$jobname'/g' \
            -e 's/^#PBS -l nodes=.*/#PBS -l nodes='$nnode':ppn='$ncpu'/g' \
            -e 's/^#PBS -l walltime=.*$/#PBS -l walltime=00:05:00/g' \
            template-job-pbs > job-$nnode.pbs
        sed -i $nheader',$d' job-$nnode.pbs

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
        echo -e "$COMMAND bmk-$nnode.namd > bmk-$nnode.log\nwait" >> job-$nnode.pbs
    done
fi
