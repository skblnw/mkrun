#! /bin/tcsh -f
########################################################
## Usage: dG_decomp.sh forward.fepout 4000 300000
##
## hard assumptions: 
##   equil steps:      4000 steps   (CHANGE IF NOT)
##   FEP steps/window: 300000 steps (CHANGE IF NOT)
##   Temperature:      300 K
##   FEPOUT needs to be written EVERY step
########################################################
##default input filename, and FEP steps, equiibration steps
set filein="alchemy.fepout"
set feptot="200000"
set fepfreq="400"

if ( $#argv ==  3 ) then
    set filein = $argv[1]
    set feptot = $argv[2]
    set fepfreq = $argv[3]
endif
echo ${filein} ${feptot} ${fepfreq}

## Find out line numbers of  each new window
set j = `awk '/#STARTING/ {print NR}'<${filein}`
#echo $j
rm -f fep_decomp.out
set id = 0 
set nlines = `expr ${feptot} / ${fepfreq} + 1`  # new NAMD2 files starts at step 1 for fepout results

foreach n($j)
  set id = `expr ${id} + 1`
  set j1 = `expr $n + 1`
  echo "window" $id $j1 $nlines
  ## fetch the correct lines in the fepout results for each window
  tail -n +${j1} ${filein} | head -${nlines} | awk '{print  ($3 - $4), ($5 - $6), $7, $9}'>tmp${id}

  ## computation of decomposition energies
  awk '{sum4 +=  $4;sum1 += (exp($1 / 0.001987 / $4));sum2 += (exp($2  / 0.001987 / $4));sum3 += (exp($3 / -0.001987 / $4));n++;} END {print ((-0.001987) * sum4 / n * log(sum1 / n)), ((-0.001987) * sum4 / n * log(sum2 / n)),((-0.001987) * sum4 / n * log(sum3 / n)),n;}'<tmp${id}>>fep_decomp.out
  rm -f tmp${id}
end

awk '{sum1 +=   $3;sum2 +=  $1; sum3 +=  $2; } END{print  "dG", sum1, "dGelec", sum2, "dGvdw", sum3;}'<fep_decomp.out
