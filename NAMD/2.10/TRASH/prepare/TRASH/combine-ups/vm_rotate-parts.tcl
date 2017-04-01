set sel [atomselect top "chain A and resid 1 to 13"]
set selcenter [atomselect top "chain A and resid 14 and name N"]
set centerX [list [$selcenter get x] [$selcenter get y] [$selcenter get z]]
set matrixX [trans center $centerX axis x 10]
$sel move $matrixX

# set sel [atomselect top "segname P3 and resid 1 to 13"]
# set selcenter [atomselect top "segname P3 and resid 14 and name N"]
# set centerX [list [$selcenter get x] [$selcenter get y] [$selcenter get z]]
# set matrixX [trans center $centerX axis y -10]
# $sel move $matrixX