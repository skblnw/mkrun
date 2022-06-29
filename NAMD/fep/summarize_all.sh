#!/bin/bash

[ $# -eq 0 ] && { echo "> Usage: $0 [anything]"; exit 1; }

for ii in `find . -name "summaryddg"`; do cat $ii; done > ddg
for ii in `find . -name "summarydg"`; do cat $ii; done > dg

# for ii in $(seq 1 9); do target=`echo GILGFVFTL | sed -e 's/./A/'$ii`; cd posi$ii; rm sum*; summarize GILGFVFTL $target; cd ..; done
