SEP_FILE=/home/kevin/Dropbox/QWD/PH-memb/smd/4/run-ss-xy/make_analysis/cal_sep/com/smd5.dat
DCD_FILE=/share/data/kevin/PH-memb/smd/4/run-ss-xy/output/smd5-1.dcd

OUTPUT_DIR=smd5
mkdir -p $OUTPUT_DIR

> $OUTPUT_DIR/FRAMES

while read line; do
    srcname=smd-20
    pdbname=us-z$line-0.pdb
    outname=us-z$line-0.restart

    nframe=`cat $SEP_FILE | awk '{print $1" "$2*10}' | egrep " $line" | head -n 1 | awk '{print $1}'`
    let nframe*=5
    while [[ $nframe -lt 1 ]]
    do
        let line+=1
        nframe=`cat $SEP_FILE | awk '{print $1" "$2*10}' | egrep " $line" | head -n 1 | awk '{print $1}'`
        let nframe*=5
        echo "Info: iterating line=$line nframe=$nframe"
    done
    
#    nstep=$(expr $nframe*100)

#    vmd -e writepdb.tcl -args $pdbname $nframe < /dev/null

    catdcd -o $OUTPUT_DIR/$outname.coor -otype namdbin -first $nframe -last $nframe -dcd $DCD_FILE
    echo $nframe >> $OUTPUT_DIR/FRAMES
#    catdcd -o $outname.vel -otype namdbin -first $nframe -last $nframe -dcd ../output/$srcname.veldcd
#    head -n 2 ../output/$srcname.xsc > header
#    grep "^$nstep " ../output/$srcname.xst > xsc
#    cat header xsc > $outname.xsc
#    rm header xsc
done < LINES

#if [ -e $srcname ]; then
#    rm -rf $srcname
#    echo "Removing folder $srcname"
#fi

#mv *.coor *.vel *.xsc $srcname
