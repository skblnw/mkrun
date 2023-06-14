#!/bin/bash

# Initialize sum and count to 0
sum=0
count=0

# Loop over all directories
for i in $(seq 1 63); do
  # Check if the directory exists
  if [ -d "t${i}" ]; then
    # Check if the log file exists
    if [ -f "t${i}/LOG_fep" ]; then
      # Get the average value from the eighth column of the grep output
      avg=$(grep Benchmark "t${i}/LOG_fep" | awk '{sum=0; count=0} {if($8!=""){sum+=$8; count++}} END {if(count>0) print sum/count; else print "NA"}')

      # If the avg is not "NA", add it to the total sum and increment the count
      if [ "$avg" != "NA" ]; then
        sum=$(echo "$sum + $avg" | bc)
        count=$((count + 1))
      fi
    fi
  fi
done

# Calculate and print the overall average, only if count is not zero
if [ $count -ne 0 ]; then
  overall_avg=$(echo "scale=2; $sum / $count" | bc)
  echo "Overall average: $overall_avg"
else
  echo "No averages were calculated."
fi
