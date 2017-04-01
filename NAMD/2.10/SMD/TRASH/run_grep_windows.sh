while read line; do

    srcname=smd-20
    pdbname=us-z$line-0.pdb
    outname=us-z$line-0.restart

    nframe=`sed '1d' ../converge/$srcname.sep | awk ' NR % 5 == 0 { print $1/100" "$2*10}' | egrep " $line" | head -n 1 | awk '{print $1}'`
    nstep=$(expr $nframe*100)

#    vmd -e writepdb.tcl -args $pdbname $nframe < /dev/null

    catdcd -o $outname.coor -otype namdbin -first $nframe -last $nframe -dcd ../output/$srcname.dcd
    catdcd -o $outname.vel -otype namdbin -first $nframe -last $nframe -dcd ../output/$srcname.veldcd
    head -n 2 ../output/$srcname.xsc > header
    grep "^$nstep " ../output/$srcname.xst > xsc
    cat header xsc > $outname.xsc
    rm header xsc

done < LINES

#if [ -e $srcname ]; then
#    rm -rf $srcname
#    echo "Removing folder $srcname"
#fi

mkdir $srcname
mv *.coor *.vel *.xsc $srcname
