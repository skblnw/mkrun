nnode=0
ncpu=8
name=500k-cpu-8omp

echo "" > bm-result-$name.dat
for ii in {2,4,6}
do
    nnode=$ii
    NS=`grep "Performance" {1..3}/500k-cpu-${nnode}node-${ncpu}omp.log | awk '{sum+=$2} END {print sum/NR}'`
    echo -e "$ii\t$NS" >> bm-result-$name.dat
done
