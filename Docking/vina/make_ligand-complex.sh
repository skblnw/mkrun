DIR_VINA=/users/PAS1326/osu9603/cpf1/vina/ligand/nat_product

for ii in $(cat list)
do
    cd $ii
    cp $DIR_VINA/template-vina-sh-complex vina-complex.sh
    sed -e 's/NAME/complex-'$ii'/g' $DIR_VINA/template-job-complex > job.pbs
    qsub job.pbs
    cd ..
done
