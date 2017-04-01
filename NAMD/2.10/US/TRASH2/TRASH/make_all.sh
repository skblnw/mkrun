#!bin/bash

winsfolder=/share/data/kevin/BAR-PH/smd-sym-ph/windows/smd-20
# Set both to 0 if it is the first time
startj=21
endj=21

count=0
colcount=0
cpcount=0

if [ $startj -eq 0 ]; then
    echo -e "As startj == 0\nInitializing..."
    mkdir output
fi

#for i in 680 695 710 725 740 755 770 785 800 815 830 845 860 875 890 905 920 935 950 965 980
for i in 616 630 640 650 660 670 680 695 710 725 740 755 765 770 785 800 815 830 845 860 875 890 905 920 935 950 965 980
#for i in 616 630 640 650 660 670
#for i in 760 765
do
    posaa=`echo "scale=1;$i / 10" | bc`

    if [ $startj -eq 0 ]; then
        if [ -e us-z$i ]; then
            echo "Opps I am deleting your folder us-z$i and making a new one"
            rm -r us-z$i
        fi
  
#        echo "Creating a folder us-z$i"
        mkdir us-z$i
  
        cd us-z$i
#        echo "Oh yeah I get into folder us-z$i"
    
        colname=colvars-z$i
        sed -e 's/POSaa/'$posaa'/g' ../template-colvars > $colname.in
#        echo "I read template-colvars and make a new one"
        colcount=$(expr $colcount + 1)

        cp $winsfolder/*$i*.* ../output
        nn=`ls $winsfolder/*$i*.* | wc -l`
        cpcount=$(expr $cpcount + $nn)

        cd ..

    else

        cd us-z$i
#        echo "Oh yeah I get into folder us-z$i"
  
        for ((j=$startj; j<=$endj; j++)); do
          k=$(expr $j - 1)
        
          fname=us-z$i-$j
          sed -e 's/INNAME/'us-z$i-$k'/g' -e 's/OUTNAME/'us-z$i-$j'/g' -e 's/POSaa/'$i'/g' ../template-namd > $fname.namd
#          echo "I read template-namd2 and make the $j-th one"
          count=$(expr $count + 1)
          gjname=us-z$i-$j
          sed -e 's/POSaa/'$i'/g' -e 's/NUM/'$j'/g' ../template-gjob > $gjname.pbs
#          echo "I read template-gjob and make the $j-th one"
          count=$(expr $count + 1)
        done
      
        cd ..

    fi
done


if [ $startj -eq 0 ]; then
    echo "$colcount colvars has been made"
    echo "$cpcount files has been copied"
fi
count=`echo "$count / 2" | bc`
echo "$count namd/pbs has been made"
