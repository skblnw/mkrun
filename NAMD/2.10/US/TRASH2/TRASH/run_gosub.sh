#!/bin/bash

# If number of windows are smaller (<12), they must be submitted one by one or sim-2 will start when sim-1 has not yet finished
for j in {12..12}
do
#  for i in 616 630 640 650 660 670 680 695 710 725 740 755 765 770 785 800 815 830 845 860 875 890 905 920 935 950 965 980
#  for i in 800 815 830 845 860 875 890 905 920 950 965 980
#  for i in 616 630 640 650 660 670 680 695 710 725 740 755 765 770 785
  for i in 616 630 640 650 660 670
  do
    foldername=./us-z$i
    gosub $foldername/us-z$i-$j.pbs
  done
done
