#!/bin/bash

if [[ "$#" -eq 0 ]]; then
    echo ">Usage: mknamd_fep_decomp alchemy.fepout [-multiple (Optional)]"
    echo ">If you don't have fep.namd"
    echo ">Usage: mknamd_fep_decomp alchemy.fepout [alchEquilSteps] [numSteps] [alchOutFreq] [-multiple (Optional)]"
    exit 0
fi

if [[ $2 == "-multiple" ]]; then
    tag_multiple=true
else
    tag_multiple=false
fi

if $tag_multiple; then 

    triallist=`find . -maxdepth 1 -name 'trial[0-9]*' | sort -V`
    for trial in $triallist; do
        [ ! -f "$trial/$1" ] && { echo -e "$trial/$1 does not exist!"; exit 1; }
        echo "Computing for $trial"
        outdir=$trial/outdecomp
        [ -d "$outdir" ] && { rm -r $outdir }
        mkdir $outdir
        grep "^#Free energy change" $trial/$1 | awk '{print NR, $9, $12, $19}' > $outdir/raw
        grep -n "^#STARTING COLLECTION" $trial/$1 | awk -F: '{print $1+1}' > $outdir/l1
        grep -n "^#Free energy change" $trial/$1 | awk -F: '{print $1}' > $outdir/l2
        paste $outdir/l1 $outdir/l2 > $outdir/l3
        awk '{print NR,$1,($2-$1)}' $outdir/l3 > $outdir/nlist
        while read -r nn ntail nhead; do 
            tail -n +$ntail $trial/$1 | head -$((nhead+1)) >> $outdir/raw$nn
            tail -n +$ntail $trial/$1 | head -$((nhead)) | awk 'NR%1==0' | awk '{print ($4 - $3), ($6 - $5), $7, $9}' >> $outdir/data$nn
        done < $outdir/nlist
    done

else

    triallist="1"
    [ ! -f "$1" ] && { echo -e "$1 does not exist!"; exit 1; }
    outdir=outdecomp
    [ -d "$outdir" ] && { rm -r $outdir }
    mkdir $outdir
    grep "^#Free energy change" $1 | awk '{print NR, $9, $12, $19}' > $outdir/raw
    grep -n "^#STARTING COLLECTION" $1 | awk -F: '{print $1+1}' > $outdir/l1
    grep -n "^#Free energy change" $1 | awk -F: '{print $1}' > $outdir/l2
    paste $outdir/l1 $outdir/l2 > $outdir/l3
    awk '{print NR,$1,($2-$1)}' $outdir/l3 > $outdir/nlist
    while read -r nn ntail nhead; do 
        tail -n +$ntail $1 | head -$((nhead+1)) >> $outdir/raw$nn
        tail -n +$ntail $1 | head -$((nhead)) | awk 'NR%1==0' | awk '{print ($4 - $3), ($6 - $5), $7, $9}' >> $outdir/data$nn
    done < $outdir/nlist

fi

rm decompose_summary.dat
for trial in $triallist; do

    if $tag_multiple; then 
        outdir=$trial/outdecomp
    else
        outdir=outdecomp
    fi

    while read -r nn v1 v2; do 
        awk '{sum_temp += $4; \
            esum_elec += (exp($1 / -0.001987 / $4)); \
            esum_vdw += (exp($2  / -0.001987 / $4)); \
            esum_total += (exp($3 / -0.001987 / $4)); \
            sum_elec += $1; sum2_elec += $1*$1; \
            sum_vdw += $2; sum2_vdw += $2*$2; \
            sum_couple += $1*$2; n++;} END \
            {printf "%3d %5.2f % 6.2f % 10.2f % 6.2f % 6.2f % 6.2f % 6.2f % 10.2f % 6.2f % 6.2f % 6.2f % 6.2f % 10.2f\n", \
                n, \
                sum_temp/n, \
                (-0.001987*sum_temp/n*log(esum_total/n)), \
                (-0.001987*sum_temp/n*log(esum_elec/n)), \
                (sum_elec/n) - (sum2_elec/n - sum_elec*sum_elec/n/n)/1.24 - (sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
                (sum_elec/n), \
                -(sum2_elec/n - sum_elec*sum_elec/n/n)/1.24, \
                -(sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
                (-0.001987*sum_temp/n*log(esum_vdw/n)), \
                (sum_vdw/n) - (sum2_vdw/n - sum_vdw*sum_vdw/n/n)/1.24 - (sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
                (sum_vdw/n), \
                -(sum2_vdw/n - sum_vdw*sum_vdw/n/n)/1.24, \
                -(sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
                -0.001987*sum_temp/n*log(esum_total/n) + 0.001987*sum_temp/n*log(esum_elec/n) + 0.001987*sum_temp/n*log(esum_vdw/n) \
                }' < $outdir/data$nn >> $outdir/decomp
        awk '{sum_temp += $4; \
            sum1 += (exp($3 / -0.001987 / $4)); \
            sum2 += (exp(2 * $3 / -0.001987 / $4)); \
            sum_total += $3; sum2_total += $3*$3; n++;} END \
            {printf "%3d %5.2f % 6.2f %6.2f % 6.2f %6.2f %6.2f %7.0f %3.2e %3.2e %3.2e\n", \
                n, \
                sum_temp/n, \
                (-0.001987*sum_temp/n*log(sum1/n)), \
                sqrt(((n*sum2/sum1/sum1-1)/0.3844)), \
                (sum_total/n), \
                sqrt(sum2_total/n - sum_total*sum_total/n/n), \
                (sum_total/n + 0.001987*sum_temp/n*log(sum1/n)), \
                (exp((sum_total/n - (-0.001987*sum_temp/n*log(sum1/n)))*n/sum_temp/0.001987)), \
                (sum1/n), \
                (sum1*sum1/n/n), \
                (sum2/n)}' < $outdir/data$nn >> $outdir/error
    done < $outdir/nlist

    paste $outdir/raw $outdir/decomp | awk '{printf "%2d % 6.2f % 6.2f % 6.2f\n", NR, $3, $7, ($3 - $7)}' > $outdir/compare_dg

    if $tag_multiple; then 
        echo "mknamd> <$trial>"
        echo -n "mknamd>    FEP gives you: "
        tail -1 $trial/$1 | awk '{print $19}'
    else
        echo -n "mknamd>    FEP gives you: "
        tail -1 $1 | awk '{print $19}'
    fi
    awk '{sum_dg += $3; \
        sum_delec += $4; \
        sum_dvdw += $9; \
        sum_couple += $14; n++;} END \
        {printf "mknamd> Decomp gives you: %.2f\n\telec   % 7.2f\n\tvdw    % 7.2f\n\tcoupl  % 7.2f\n", sum_dg, sum_delec, sum_dvdw, sum_couple}' < $outdir/decomp
    awk -v "prefix=${trial#./trial}" '{sum_dg += $3; \
        sum_delec += $4; \
        sum_dvdw += $9; \
        sum_couple += $14; n++;} END \
        {print prefix,sum_dg,sum_delec,sum_dvdw,sum_couple}' < $outdir/decomp >> decompose_summary.dat

done

if $tag_multiple; then 
    echo -ne "mknamd> Average over all trials\nmknamd> dG dElec dVdw dCoupl\nmknamd> "
    awk '{sum_dg += $2; \
        sum_delec += $3; \
        sum_dvdw += $4; \
        sum_couple += $5; n++;} END \
        {printf "%.2f %.2f %.2f %.2f", sum_dg/n,sum_delec/n,sum_dvdw/n,sum_couple/n}' < decompose_summary.dat
    echo ""
fi
