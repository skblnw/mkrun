#!/bin/bash

DIR=output
echo -e "Target directory is [$DIR]"
filename="NPT"
echo -e "File name pattern is [$filename]"
mkdir -p output-dcd output-restart
FIRST=`ls output/*.dcd | sort -V | head -1 | awk -F '-' '{print $2}' | awk -F '.' '{print $1}'`
LAST=`ls output/*.dcd | sort -V | tail -1 | awk -F '-' '{print $2}' | awk -F '.' '{print $1}'`
echo -e "The smallest/first (version) number of DCD files is $FIRST"
echo -e "The largest/latest (version) number of DCD files is $LAST"

# Prompted to ask if user want to continue
read -p "Start moving?" -n 1 -r
echo 
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Stopped"
    exit 1
fi

mv $DIR/$filename-$LAST.restart.* output-restart
eval mv $DIR/$filename-{$FIRST..$LAST}.dcd output-dcd
