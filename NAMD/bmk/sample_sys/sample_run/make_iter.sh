for nn in {1..2}
do
    mkdir -p $nn
    mkdir -p output$nn
    sed -e 's/NN/'$nn'/g' make_prepare.sh > tmp.sh
    sh tmp.sh
    sh make_submit.sh
done
