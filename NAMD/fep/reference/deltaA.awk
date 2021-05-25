#!/bin/awk -f
################################################################################
# Extract free energy differences from alchemical FEP in NAMD 2.5b3 and later  #
#                                                                              #
# Example : deltaA.awk *.fepout > fep.dat                                      #
#                                                                              #
# Computes free energy as a function of lambda                                 #
# Computes average electrostatic energy difference for component analysis in TI#
#                                                                              #
# Jerome Henin <jerome.henin@uhp-nancy.fr>                                     #
################################################################################

BEGIN {
    A = 0
    Eel = 0
    FLAG_L2 = 0
    start = 1
    printf "# lambda         A          dA       sum<dEel>      <dEel>\n"
}

$1=="#Free" {
    L = $8
    L2 = $9
    if (FLAG_L2 && prev_L2 != L) {
        printf "# Non-contiguous windows: %f --> %f\n", prev_L2, L
        A = 0
        Eel = 0
        FLAG_L2 = 0
        start = 1
    }
    if ( start ) {
        start = 0;
        printf "%-8s  %9.3f    %9.3f    %9.3f    %9.3f\n", L, A, 0, 0, 0
    }
    prev_L2 = L2
    FLAG_L2 = 1
    dA = $12
    A += dA
    if (n>0)  {
        dEel /= n
        Eel += dEel
    }
    n = 0
    printf "%-8s  %9.3f    %9.3f    %9.3f    %9.3f\n", L2, A, dA, Eel, dEel
}

$1=="FepEnergy:" {
    dEel += $4 - $3
    n++
}

