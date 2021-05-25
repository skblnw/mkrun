#!/bin/bash

if [[ "$#" -ne 1 ]]; then
    echo ">Usage: mknamd_fep_grep alchemy.fepout"
    exit 0
fi

grep "^#Free energy change" $1 | awk '{print NR" "$12" "$19}' > fepout

cat > pyplot <<EOF
import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('fepout', comments=['#','@'])
plt.figure(1)
ax = plt.gca()
ax.bar(data[:,0], data[:,1])
plt.figure(2)
ax = plt.gca()
ax.plot(data[:,0], data[:,2], linewidth=2, marker="o")
plt.show()
EOF

python pyplot
rm -f pyplot fepout