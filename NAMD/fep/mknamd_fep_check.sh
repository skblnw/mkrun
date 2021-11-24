#!/bin/bash

if [[ "$#" -eq 0 ]]; then
    echo ">Usage: mknamd_fep_grep alchemy.fepout <-plot> <-range start end>"
    exit 0
fi

if [[ "$2" == "-range" ]]; then 

    for ii in $(seq $3 $4); do 
        [ ! -f "t$ii/$1" ] && { echo -e "t$ii/$1 does not exist!"; exit; }
        grep "^#Free energy change" t$ii/$1 | awk '{print NR" "$9" "$12" "$19}' > $ii
    done

else

    [ ! -f "$1" ] && { echo -e "$1 does not exist!"; exit; }
    grep "^#Free energy change" $1 | awk '{print NR" "$9" "$12" "$19}' > fepout
    
fi

if [[ "$2" == "-plot" ]]; then 

    cat > pyplot <<EOF
import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('fepout', comments=['#','@'])
plt.figure(1)
ax = plt.gca()
ax.bar(data[:,0], data[:,2])
plt.figure(2)
ax = plt.gca()
ax.plot(data[:,0], data[:,-1], linewidth=2, marker="o")
plt.show()
EOF

    python pyplot
    rm -f pyplot fepout
    
fi
