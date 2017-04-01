# mol new protein_nolinker_coil-rotated_centered_autopsf.psf
# mol addfile protein_nolinker_coil-rotated_centered_autopsf.pdb

# mol new membrane_duplicated.psf
# mol addfile membrane_final_unwrap_centered.pdb

set sel1 [atomselect 0 all]
set sel2 [atomselect 1 all]

package require topotools

set mol [::TopoTools::selections2mol "$sel1 $sel2"]
set selall [atomselect $mol all]

animate write psf system.psf $mol
animate write pdb system.pdb $mol