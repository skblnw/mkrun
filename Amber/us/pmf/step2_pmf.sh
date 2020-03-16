#!/bin/bash

wham 32 54.5 48 0.00000001 303 0 metadata out.pmf
sed '1d' out.pmf | awk '{print $1,"",$2}' > plot_free_energy.dat
echo "Now xmgrace plot_free_energy.dat"
