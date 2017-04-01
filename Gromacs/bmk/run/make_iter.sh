for nn in {1..1}
do
    mkdir -p $nn
    sed -e 's/NUM_ITER/'$nn'/g' make_prepare.sh > tmp.sh
    sh tmp.sh
    sh make_submit.sh
done
