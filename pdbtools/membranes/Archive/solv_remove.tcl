package require solvate
solvate 5kxi_popc.psf 5kxi_popc.pdb -o 5kxi_popc_water_TEMP -b 1.5 -minmax {{-38 -38 -39} {39 39 50}}
set all [atomselect top all]
$all set beta 0
set seltext "segid WT1 to WT99 and same residue as abs(z) <25"
set sel [atomselect top $seltext]
$sel set beta 1
set badwater [atomselect top "name OH2 and beta > 0"]
set seglist [$badwater get segid]
set reslist [$badwater get resid]
mol delete all
package require psfgen
resetpsf
readpsf 5kxi_popc_water_TEMP.psf
coordpdb 5kxi_popc_water_TEMP.pdb
foreach segid $seglist resid $reslist {
delatom $segid $resid
}
writepdb 5kxi_popcw.pdb
writepsf 5kxi_popcw.psf
file delete kcsa_popc_water_TEMP.psf
file delete kcsa_popc_water_TEMP.pdb

