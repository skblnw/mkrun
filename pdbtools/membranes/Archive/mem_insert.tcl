set POPC "resname POPC"
set all [atomselect top all]
$all set beta 0
set seltext1 "$POPC and same residue as (within 0.6 of protein)"
set sel1 [atomselect top $seltext1]
$sel1 set beta 1
set badlipid [atomselect top "name P and beta > 0"]
$badlipid get resid
set sel [atomselect top "beta<1 "]
$sel writepdb prot_memb_hole.pdb
