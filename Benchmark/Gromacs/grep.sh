#!/bin/bash

filename="bm-result"
rm -f $filename
touch $filename
for ii in {1,2,3,4,8,10}
do
NS=`grep -i "Performance" $ii.log | awk '{print $2}'`
echo -e "$ii\t$NS" >> $filename
done
