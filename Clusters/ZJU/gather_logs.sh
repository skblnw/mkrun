#!/bin/bash

# Specify the start and end date
start_date="20230507"
end_date="20230522"

# Output CSV header
echo "Node,Date,Time,GPU,Temperature,Utilization" > all_logs.csv

# Specify the directory containing the logs
dir="gpu_temps_logs"

# Process log files sorted by node number
for file in $(ls $dir/node*_*.log | sort -V); do
    # Extract node from filename
    node=$(basename $file | cut -d'_' -f1)

    # Extract date from filename and remove file extension
    date=$(basename $file | cut -d'_' -f2 | cut -d'.' -f1)
    
    # Check if date falls within the desired range
    if [[ "$date" < "$start_date" ]] || [[ "$date" > "$end_date" ]]; then
        continue
    fi

    # Read and process the entire log file at once
    awk -v node="$node" -v date="$date" '{
        split($1,time,":");
        print node","date","time[1]","$2","$3","$4
    }' $file >> all_logs.csv
done
