#!/bin/bash

rm -f 0mini *eq 5tgtus
for ii in $(seq 1 32)
do
    echo pmemd.cuda -O -i win$ii/0mini.mdin.$ii -p win$ii/ionized.parm7.$ii -c win$ii/ionized.rst7.$ii -ref win$ii/ionized.rst7.$ii \
    -x win$ii/0mini.x.$ii -o win$ii/0mini.out.$ii -inf win$ii/0mini.info.$ii -r win$ii/0mini.r.$ii \& >> 0mini
    for neq in $(seq 1 4)
    do
    if [ $neq -eq 1 ]; then
        neq_1=0mini
    else
        jj=$(($neq - 1))
        neq_1=${jj}eq
    fi
    echo -O -i win$ii/${neq}eq.mdin.$ii -p win$ii/ionized.parm7.$ii -c win$ii/${neq_1}.r.$ii -ref win$ii/ionized.rst7.$ii \
    -x win$ii/${neq}eq.x.$ii -o win$ii/${neq}eq.out.$ii -inf win$ii/${neq}eq.info.$ii -r win$ii/${neq}eq.r.$ii >> ${neq}eq
    done
    echo -O -i win$ii/5tgtus.mdin.$ii -p win$ii/ionized.parm7.$ii -c win$ii/${neq}eq.r.$ii -ref win$ii/ionized.rst7.$ii \
    -x win$ii/5tgtus.x.$ii -o win$ii/5tgtus.out.$ii -inf win$ii/5tgtus.info.$ii -r win$ii/5tgtus.r.$ii >> 5tgtus
done

cat > 0mini.mdin << EOF
Minimization input file in explicit solvent
&cntrl
 imin=1,maxcyc=5000,ncyc=2500,
 cut=12.0,fswitch=10.0,
 ntpr=100,ntxo=2,
 ntr=1,restraintmask='!:WAT & !@H=,K+,Na+,Cl-',restraint_wt=10,
 watnam='WAT',owtnm='O',     
/
&ewald vdwmeth = 0 /
&wt type='END' /
LISTIN=POUT
LISTOUT=POUT
&end
EOF

cat > 1eq.mdin << EOF
A NVT simulation for common production-level simulations
&cntrl
 imin=0,irest=0,ntx=1,
 ntt=3,gamma_ln=1.0,tempi=303.15,temp0=303.15,
 cut=12.0,fswitch=10.0,
 nstlim=125000,dt=0.001,
 ntc=2,ntf=2,
 ntpr=1000,ntwx=12500,ntwr=12500,

 ntxo=2,ioutfm=1,iwrap=0,
 ntr=1,restraintmask='!:WAT & !@H=,K+,Na+,Cl-',restraint_wt=10,
 watnam='WAT',owtnm='O',
 infe=1,
/ 
&ewald vdwmeth = 0 /
&wt type='END' /
&pmd output_file='winNN/1eq.pmd.NN',output_freq=50,cv_file='winNN/cv.in' /
LISTIN=POUT
LISTOUT=POUT
&end
EOF

cat > 2eq.mdin << EOF
A NVT simulation for common production-level simulations
&cntrl
 imin=0,irest=1,ntx=5,
 ntt=3,gamma_ln=1.0,temp0=303.15,
 cut=12.0,fswitch=10.0,
 nstlim=125000,dt=0.001,
 ntc=2,ntf=2,
 ntpr=1000,ntwx=12500,ntwr=12500,

 ntxo=2,ioutfm=1,iwrap=0,
 ntr=1,restraintmask='!:WAT & !@H=,K+,Na+,Cl-',restraint_wt=10,
 watnam='WAT',owtnm='O',    
/ 
&ewald vdwmeth = 0 /
&wt type='END' /
&pmd output_file='winNN/2eq.pmd.NN',output_freq=50,cv_file='winNN/cv.in' /
LISTIN=POUT
LISTOUT=POUT
&end
EOF

cat > 3eq.mdin << EOF
A NPT simulation for common production-level simulations
&cntrl
 imin=0,irest=1,ntx=5,
 ntt=3,gamma_ln=1.0,temp0=303.15,
 cut=12.0,fswitch=10.0,
 nstlim=125000,dt=0.001,
 ntc=2,ntf=2,
 ntpr=1000,ntwx=12500,ntwr=12500,

 ntxo=2,ioutfm=1,iwrap=0,
 barostat=2,ntp=1,pres0=1.0,
 ntr=1,restraintmask='!:WAT & !@H=,K+,Na+,Cl-',restraint_wt=10,
 watnam='WAT',owtnm='O',  
/ 
&ewald vdwmeth = 0 /
&wt type='END' /
&pmd output_file='winNN/3eq.pmd.NN',output_freq=50,cv_file='winNN/cv.in' /
LISTIN=POUT
LISTOUT=POUT
&end
EOF

cat > 4eq.mdin << EOF
A NPT simulation for common production-level simulations
&cntrl
 imin=0,irest=1,ntx=5,
 ntt=3,gamma_ln=1.0,temp0=303.15,
 cut=12.0,fswitch=10.0,
 nstlim=250000,dt=0.001,
 ntc=2,ntf=2,
 ntpr=1000,ntwx=12500,ntwr=12500,

 ntxo=2,ioutfm=1,iwrap=0,
 barostat=2,ntp=1,pres0=1.0,
 ntr=1,restraintmask='!:WAT & !@H=,K+,Na+,Cl-',restraint_wt=1,
 watnam='WAT',owtnm='O',  
/ 
&ewald vdwmeth = 0 /
&wt type='END' /
&pmd output_file='winNN/4eq.pmd.NN',output_freq=50,cv_file='winNN/cv.in' /
LISTIN=POUT
LISTOUT=POUT
&end
EOF

cat > 5tgtus.mdin << EOF
A NPT simulation for common production-level simulations
&cntrl
 imin=0,irest=1,ntx=5,
 ntt=3,gamma_ln=1.0,temp0=303.15,
 cut=12.0,fswitch=10.0,
 nstlim=250000,dt=0.001,
 ntc=2,ntf=2,
 ntpr=1000,ntwx=12500,ntwr=12500,

 ntxo=2,ioutfm=1,iwrap=0,
 barostat=2,ntp=1,pres0=1.0,
 ntr=0,
 watnam='WAT',owtnm='O',  
/ 
&ewald vdwmeth = 0 /
&wt type='END' /
&pmd output_file='winNN/5tgtus.pmd.NN',output_freq=50,cv_file='winNN/cv.in' /
LISTIN=POUT
LISTOUT=POUT
&end
EOF

for ii in $(seq 1 32)
do
    for jj in 0mini 1eq 2eq 3eq 4eq 5tgtus; do
        sed -e "s/NN/$ii/g" $jj.mdin > win$ii/$jj.mdin.$ii
    done
done
