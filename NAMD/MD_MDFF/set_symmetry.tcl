set selall [atomselect top all]
$selall set occupancy 0
$selall set beta 0
set sel [atomselect top "name CA"]
$sel set occupancy 1
set sel1 [atomselect top "(segname AP1 or segname BP1) and name CA"]
$sel1 set beta 1
set sel2 [atomselect top "(segname CP1 or segname DP1) and name CA"]
$sel2 set beta 2
set sel3 [atomselect top "(segname EP1 or segname FP1) and name CA"]
$sel3 set beta 3
set sel4 [atomselect top "(segname GP1 or segname HP1) and name CA"]
$sel4 set beta 4
set sel5 [atomselect top "(segname IP1 or segname JP1) and name CA"]
$sel5 set beta 5
set sel6 [atomselect top "(segname KP1 or segname LP1) and name CA"]
$sel6 set beta 6
set sel7 [atomselect top "(segname MP1 or segname NP1) and name CA"]
$sel7 set beta 7
set sel8 [atomselect top "(segname OP1 or segname PP1) and name CA"]
$sel8 set beta 8
set sel9 [atomselect top "(segname QP1  or segname RP1) and name CA"]
$sel9 set beta 9
$selall writepdb helix-symmetry.pdb
