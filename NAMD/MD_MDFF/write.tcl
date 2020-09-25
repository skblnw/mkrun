set PREFIX ysingle
mol new ${PREFIX}_charmmgui.psf
mol addfile ${PREFIX}_charmmgui.pdb
foreach ii {1 3 5} {mol addfile $ii.coor}
set sel [atomselect top all]
foreach ii {1 2 3} jj {1 3 5} {$sel frame $ii; $sel writepdb $jj.pdb}
quit
