set filename ../combine/system
package require solvate
solvate $filename.psf $filename.pdb -b 2.4 -minmax {{-146 -72 -55} {146 72 115}} -o solvate_TEMP
mol new solvate_TEMP.psf
mol addfile solvate_TEMP.pdb
set all [atomselect top all]
$all set beta 0
set seltext "segid WT1 to WT99 and same residue as abs(z) < 20"
set sel [atomselect top $seltext]
$sel set beta 1
set badwater [atomselect top "name OH2 and beta > 0"]
set seglist [$badwater get segid]
set reslist [$badwater get resid]
mol delete all
package require psfgen
resetpsf
readpsf solvate_TEMP.psf
coordpdb solvate_TEMP.pdb
foreach segid $seglist resid $reslist {
delatom $segid $resid
}
writepdb solvated.pdb
writepsf solvated.psf
#file delete solvated_TEMP.psf
#file delete solvated_TEMP.pdb

package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.18 -cation POT -o ionized
quit
