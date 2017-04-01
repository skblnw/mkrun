filename="bm-result"
echo "" > $filename.dat
for ii in {1,2,3,4,8,12}
do
DAY=`grep -i "benchmark time" {1..10}/bmk-$ii.log | awk '{sum+=$8} END {print sum/NR}'`
NS=`echo "scale=2; 1.0/$DAY" | bc`
echo -e "$ii\t$NS" >> $filename.dat
done
