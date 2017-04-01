RUNFOLDER=.
OUTPUTFOLDER=../output
for ii in {0..0}
do
    count=0
    for dd in 860 875 890 905 920 935 950 965 980 990 1000
    do
        list="coor xsc"
        for extension in $list; do
          if [[ ! -s $OUTPUTFOLDER/us-z$dd-$ii.restart.$extension ]]; then
              echo "us-z$dd-$ii.restart.$extension.old does not exist"
              count=$(expr $count + 1)
          fi
        done
    done
    echo "Missing:$count"
done
