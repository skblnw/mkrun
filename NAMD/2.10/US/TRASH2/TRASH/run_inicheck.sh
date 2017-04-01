#!/bin/bash

count=0

for i in 680 695 710 725 740 755 770 785 800 815 830 845 860 875 890 905 920 935 950 965 980; do
  RUNFOLDER=./us-z$i
  OUTPUTFOLDER=./output

  list="in"
  for extension in $list; do
    if [ ! -s $RUNFOLDER/colvars-z$i.$extension ]; then
        echo "NOT EXIST: $i, $extension"
        count=$(expr $count + 1)
    fi
  done

  list="coor vel xsc"
  for extension in $list; do
    if [[ ! -s $OUTPUTFOLDER/us-z$i-0.restart.$extension ]]; then
        echo "NOT EXIST: $i, $extension"
        count=$(expr $count + 1)
    fi
  done

  list="namd pbs"
  for extension in $list; do
    if [ ! -s $RUNFOLDER/us-z$i-1.$extension ]; then
        echo "NOT EXIST: $i, $extension"
        count=$(expr $count + 1)
    fi
  done
done

list="main.pdb ref.pdb ionized.pdb ionized.psf"
for name in $list; do
    if [[ ! -s $name ]]; then
        echo "NOT EXIST: $name"
        count=$(expr $count + 1)
    fi
done

if [ $count -eq 0 ]; then
    echo "Suitable for the run"
else
    echo "$count files not exist"
fi
