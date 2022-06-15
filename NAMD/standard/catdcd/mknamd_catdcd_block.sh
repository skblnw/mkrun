#!/bin/bash

if [ $# -eq 0 ]; then
    echo "make_XXX </path/to/trajectory/NPT-{1..N}.dcd>"
    exit 1
fi

end=`catdcd -num $@ | grep "Total frames" | awk '{print $3}'`

echo -e "There are in total $end frames."
read -p $'Do you want an offset? [0]\n> ' NO
if ! [[ $NO =~ ^-?[0-9]+$ ]]; then
    echo "Must be an integer"
    exit 1
fi
end=$(($end-$NO))
echo -e "There are in total $end frames."
read -p $'How many ps per frame? [1000]\n> ' PSPF
if ! [[ $PSPF =~ ^-?[0-9]+$ ]]; then
    echo "Must be an integer"
    exit 1
fi

read -p $'How many blocks? [10]\n> ' NB
if ! [[ $NB =~ ^-?[0-9]+$ ]]; then
    echo "Must be an integer"
    exit 1
fi

read -p $'How many frame for one block? [100]\nBe careful to the interval of the input DCD!\n> ' NF
if ! [[ $NF =~ ^-?[0-9]+$ ]]; then
    echo "Must be an integer\n"
    exit 1
fi

NNS=`echo "$PSPF * $NF / 1000" | bc `

mkdir -p ${NNS}ns
for (( ii=1; ii<=$NB; ii++ ))
do
    start=$(($end-$NF+1))
    echo -e "start: $start: end: $end"
    catdcd -o ${NNS}ns/block$ii.dcd -first $start -last $end $@
    wait
    end=$(($start-1))
done
