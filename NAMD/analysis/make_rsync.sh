#!/bin/bash

CUR_DIR=system-monomer-ups
echo -e "You wanna move directory [ $CUR_DIR ]"
TAR_DIR=/media/P3/kevin/Ups/system-monomer-ups
echo -e "Target directory is [ $TAR_DIR ]"
unset EXCLUDE_DIR
EXCLUDE_DIR=output/

if [ -z ${EXCLUDE_DIR+x} ]; then
    echo -e "Exclude directory [ *NULL* ]"
else
    echo -e "Exclude directory [ $EXCLUDE_DIR ]"
fi

# Prompted to ask if user want to continue
read -p "Start moving?" -n 1 -r
echo 
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Stoppped"
    exit 1
fi

if [ -z ${EXCLUDE_DIR+x} ]; then
    nohup rsync -avh ./$CUR_DIR $TAR_DIR | tee LOG_rsync-1.log &
else
    nohup rsync -avh --exclude=$EXCLUDE_DIR ./$CUR_DIR $TAR_DIR | tee LOG_rsync-1.log &
fi
