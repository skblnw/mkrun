#!/bin/bash
GMX="gmx"
MDRUN="gmx mdrun"
#MDRUN="gmx mdrun -ntmpi 1 -ntomp 12"

for ii in $(seq 1 1); do
    prefix1=step4_eq_npt
    prefix2=t${ii}
    rm -f $prefix2.tpr
    $GMX grompp -f step5_md.mdp -o output/$prefix2.tpr -c output/$prefix1.gro -r pdb2gmx/ionized.pdb -n pdb2gmx/index.ndx -p pdb2gmx/topol.top
    
    $MDRUN -v -deffnm output/$prefix2 

#    jj=$(expr $ii - 1)
#    export CUDA_VISIBLE_DEVICES=$jj
#    $MDRUN -v -deffnm output/$prefix2 > LOG$jj 2>&1 &
done
