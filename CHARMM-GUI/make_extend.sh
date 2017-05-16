for ii in {2..3}
do
    jj=$(expr $ii - 1)
    sed \
        -e 's/^set inputname.*$/set inputname         step7.'$jj'_production;/g' \
        -e 's/^outputName.*$/outputName         step7.'$ii'_production;/g' \
        template-production > step7.${ii}_production.inp
done
