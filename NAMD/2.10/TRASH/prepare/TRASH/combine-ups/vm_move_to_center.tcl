set selp [atomselect 0 all]
set selmem [atomselect 1 all]
set A [measure center $selmem]
$selmem moveby [vecinvert $A]
set B [vecsub [measure center $selmem] [measure center $selp]]
$selp moveby $B