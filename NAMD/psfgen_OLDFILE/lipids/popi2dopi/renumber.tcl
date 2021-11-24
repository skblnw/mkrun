mol new test.psf
mol addfile test.pdb

set sel [atomselect top all] 
$sel set resid [$sel get residue] 
$sel writepdb renumbered.pdb
