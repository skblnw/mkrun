DIR_VINA=/users/PAS1326/osu9603/cpf1/vina/ligand/nat_product_refined/

for rec in $(seq 2 4)
do
    for ii in $(seq 1 5)
    do
        cd $ii
        sed -e 's/NN/'$rec'/g' $DIR_VINA/template-vina-apo > vina-apo-rec$rec.sh
        sed -e 's/NN/'$rec'/g' \
            -e 's/NAME/apo-rec'$rec'-'$ii'/g' $DIR_VINA/template-job > job-rec$rec.pbs
        qsub job-rec$rec.pbs
        cd ..
    done
done
