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
##default input filename, equiibration steps, and FEP steps 
set filein="binding.fepout"
set fepeq="4000"
set feptot="300000"

if ( $#argv ==  3 ) then
    set filein = $argv[1]
    set fepeq = $argv[2]
    set feptot = $argv[3]
endif
echo ${filein} ${fepeq} ${feptot}

## Find out line numbers of  each new window
set j=`awk '/#NEW/ {print NR}'<${filein}`
echo $j
rm -f fep_decomp.out
set id = 0 
set idx = `expr ${fepeq} + 3`      # the 3 head lines in the beginning
#set idy = `expr ${feptot} - ${fepeq} +1`  # some early files starts at step 0 for fepout results
set idy = `expr ${feptot} - ${fepeq}`  # new NAMD2 files starts at step 1 for fepout results

foreach n($j)
  set id = `expr ${id} + 1`
  rm -f test_fep${id}
  set j1=`expr $n + ${idx}`
  echo "window" $id $j1 $idy
  ## fetch the correct lines in the fepout results for each window
  tail +${j1} ${filein}|head -${idy}| awk '{print  ($3 - $4), ($5 - $6), $7, $9}'>test_fep${id}

  ## computation of decomposition energies
  awk '{sum3 +=  $4;sum += (exp($3 / -0.001987 / $4));sum1 += (exp($1 / 0.001987 / $4));sum2 += (exp($2  / 0.001987 / $4));n++;} END {print ((-0.001987) * sum3 / n * log(sum1 / n)), ((-0.001987) * sum3 / n * log(sum2 / n)),((-0.001987) * sum3 / n * log(sum / n)),n;}'<test_fep${id}>>fep_decomp.out
  rm -f test_fep${id}
end

awk '{sum1 +=   $3;sum2 +=  $1; sum3 +=  $2; } END{print  "dG", sum1, "dGelec", sum2, "dGvdw", sum3;}'<fep_decomp.out
