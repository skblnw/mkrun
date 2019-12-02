for ii in $(seq 104 110); do
    sed -e 's/XXX/'$ii'/g' run-umbrella.sh > frame-${ii}_run-umbrella.sh
done