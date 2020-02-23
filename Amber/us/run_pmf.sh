#!/bin/bash

sed '1d' out.pmf | awk '{print $1,"",$2}' > plot_free_energy.dat
