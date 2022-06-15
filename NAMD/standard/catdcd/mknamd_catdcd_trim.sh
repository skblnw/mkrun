#!/bin/bash

if [ $# -eq 0 ]; then
    echo "make_XXX <prefix>(NPT) <first> <last> <stride>(5) <final ps>(10)"
    exit 1
elif [ $# -ne 5 ]; then
    echo "Wrong argument number"
    echo "make_XXX <prefix>(NPT) <first> <last> <stride>(5) <final ps>(10)"
    exit 1
fi

echo "#!/bin/bash

prefix=$1
first=$2
first_1=\`expr \$first - 1\`
last=$3
last_1=\`expr \$last - 1\`
stride=$4
ps=$5

for ii in 1
do
    unset filename
    for ((nn=\$first; nn <= \$last; nn++)); do
        if [ ! -f \$prefix-\$nn.dcd ]; then
            echo "\$prefix-\$nn.dcd NOT found!"
            exit 1
        fi

        filename=\"\$filename \$prefix-\$nn.dcd\"
    done

    if [ \$first -eq 1 ]; then
        catdcd -o \$prefix-\$last-pf\${ps}ps.dcd -stride \$stride \$filename >> LOG_trim_\$ii.log
    else
        if [ ! -f \$prefix-\${first_1}-pf\${ps}ps.dcd ]; then
            echo "Previous DCD NOT found!"
            exit 1
        fi
        catdcd -o \$prefix-tmp.dcd -stride \$stride \$filename >> LOG_trim_\$ii.log
        catdcd -o \$prefix-\$last-pf\${ps}ps.dcd \$prefix-\${first_1}-pf\${ps}ps.dcd \$prefix-tmp.dcd >> LOG_trim_\$ii.log
    fi

    rm -fv \$filename >> LOG_trim_\$ii.log
done

" > catdcd_trim.sh
