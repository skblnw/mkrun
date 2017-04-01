package require ssrestraints
ssrestraints -psf ../ionized.psf -pdb ../ionized.pdb -o extrabonds.txt -hbonds

mol new ../ionized.psf
mol addfile ../ionized.pdb
package require cispeptide
cispeptide restrain -o extrabonds-cispeptide.txt
package require chirality
chirality restrain -o extrabonds-chirality.txt

quit
