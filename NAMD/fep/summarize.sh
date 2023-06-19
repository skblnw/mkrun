#!/bin/bash

OUTPUT=summary
[ $# -ne 2 ] && { echo "> Usage: $0 [original] [target]"; exit 1; }
original="$1"
target="$2"

[ ! -f "bound/decompose_summary.dat" ] && { echo -e "bound/decompose_summary.dat does not exist!"; exit 1; }
[ ! -f "free/decompose_summary.dat" ] && { echo -e "free/decompose_summary.dat does not exist!"; exit 1; }

rm -f ${OUTPUT}dg ${OUTPUT}ddg

echo -n "$original,$target," > ${OUTPUT}ddg
awk '{sum_dg += $2; sum2_dg += $2*$2;\
      sum_delec += $3; sum2_delec += $3*$3;\
      sum_dvdw += $4; sum2_dvdw += $4*$4;\
      sum_couple += $5; sum2_couple += $5*$5; n++;} END \
      {printf "%.2f %.2f %.2f %.2f %.2f %.2f \n", sum_dg/n,sqrt((sum2_dg*n-sum_dg*sum_dg)/(n*n*(n-1))),sum_delec/n,sqrt((sum2_delec*n-sum_delec*sum_delec)/(n*n*(n-1))),sum_dvdw/n,sqrt((sum2_dvdw*n-sum_dvdw*sum_dvdw)/(n*n*(n-1)))}' < bound/decompose_summary.dat > tmp

awk '{sum_dg += $2; sum2_dg += $2*$2;\
      sum_delec += $3; sum2_delec += $3*$3;\
      sum_dvdw += $4; sum2_dvdw += $4*$4;\
      sum_couple += $5; sum2_couple += $5*$5; n++;} END \
      {printf "%.2f %.2f %.2f %.2f %.2f %.2f \n", -sum_dg/n,sqrt((sum2_dg*n-sum_dg*sum_dg)/(n*n*(n-1))),-sum_delec/n,sqrt((sum2_delec*n-sum_delec*sum_delec)/(n*n*(n-1))),-sum_dvdw/n,sqrt((sum2_dvdw*n-sum_dvdw*sum_dvdw)/(n*n*(n-1)))}' < free/decompose_summary.dat >> tmp

awk '{ddg += $1;
      ddelec += $3;
      ddvdw += $5;
      err2g += $2*$2;
      err2elec += $4*$4;
      err2vdw += $6*$6; n++;} END \
      {printf "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f \n", ddg,sqrt(err2g),ddelec,sqrt(err2elec),ddvdw,sqrt(err2vdw)}' < tmp >> ${OUTPUT}ddg

for prefix in free bound; do
  while read -r line; do
    echo -n "$original,$target,$prefix," >> ${OUTPUT}dg
    echo $line | awk '{print $2","$3","$4","$5",1,"1}' >> ${OUTPUT}dg
  done < $prefix/decompose_summary.dat
  if [ -f "$prefix/Archive/decompose_summary.dat" ]; then
    while read -r line; do
      echo -n "$original,$target,$prefix," >> ${OUTPUT}dg
      echo $line | awk '{print $2","$3","$4","$5",1,"0}' >> ${OUTPUT}dg
    done < $prefix/Archive/decompose_summary.dat
  fi
done

rm tmp
