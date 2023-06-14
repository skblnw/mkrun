#!/bin/bash

# Find all the target files
target_files=$(find . -wholename "*dir0/md.log")

# Use awk to get the number from the second column, sum them up and count the number of lines
result=$(awk '/Performance/ {sum+=$2; count++} END {if(count>0) print sum/count; else print "No valid data found."}' $target_files)

echo "The average performance value is: $result"
