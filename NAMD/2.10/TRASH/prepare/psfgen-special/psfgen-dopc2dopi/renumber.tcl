mol new membrane_psfgen.psf
mol addfile membrane_psfgen.pdb

set sel [atomselect top all] 
$sel set chain L
$sel set resid [$sel get residue]
$sel writepdb renumbered.pdb
$sel writepsf renumbered.psf