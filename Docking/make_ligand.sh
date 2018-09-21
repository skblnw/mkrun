DIR_VINA=/users/PAS1326/osu9603/cpf1/vina/ligand/pbc_pdbqt

for ii in $(seq 1 2)
do
    cd $ii
    cp $DIR_VINA/template-vina-sh vina-apo.sh
    sed -e 's/NAME/apo-'$ii'/g' $DIR_VINA/template-job > job.pbs
    qsub job.pbs
    cd ..
done
