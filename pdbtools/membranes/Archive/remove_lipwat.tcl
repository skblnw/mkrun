mol load psf 5kxi_popc_raw.psf pdb 5kxi_popc_raw.pdb
set POPC "resname POPC"
set all [atomselect top all]
$all set beta 0
set seltext1 "$POPC and same residue as (name P1 and z>0 and abs(x)<15 and abs(y)<15)"
set seltext2 "$POPC and same residue as (name P1 and z<0 and abs(x)<10 and abs(y)<10)"
set seltext3 "$POPC and same residue as (within 0.6 of protein)"
set sel1 [atomselect top $seltext1]
set sel2 [atomselect top $seltext2]
set sel3 [atomselect top $seltext3]
$sel1 set beta 1
$sel2 set beta 1
$sel3 set beta 1
set badlipid [atomselect top "name P1 and beta >0"]
set seglistlipid [$badlipid get segid]
set reslistlipid [$badlipid get resid]
set seltext4 "(water and not segname WCA WCB WCC WCD WF SOLV) and same residue as within 3 of ((same residue as (name P1 and beta>0)) or protein)"
set seltext5 "segname SOLV and same residue as within 3 of lipids"
set sel4 [atomselect top $seltext4]
set sel5 [atomselect top $seltext5]
$sel4 set beta 1
$sel5 set beta 1
set badwater [atomselect top "name OH2 and beta >0"]
set seglistwater [$badwater get segid]
set reslistwater [$badwater get resid]
mol delete all
package require psfgen
resetpsf
readpsf 5kxi_popc_raw.psf
coordpdb 5kxi_popc_raw.pdb
foreach segid $seglistlipid resid $reslistlipid {
delatom $segid $resid
}
foreach segid $seglistwater resid $reslistwater {
delatom $segid $resid
}
writepsf 5kxi_popc.psf
writepdb 5kxi_popc.pdb

