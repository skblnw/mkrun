for win in 616 630 640 650 660 670 680 695 710 725 740 755 765 770 785 800 815 830 845 860 875 890 905 920 935 950 965 980
do
    catdcd -o /mnt/flash/kevin/data2/us-sym/output/us-z$win-1-10.dcd -stride 5 output/us-z$win-{1..10}.dcd >& LOG_catdcd.log
done
