set filename ../prot_psfgen
package require solvate
#solvate $filename.psf $filename.pdb -b 2.4 -minmax {{-77 -76 -35} {77 76 130}} -o solvated
solvate $filename.psf $filename.pdb -b 2.4 -t 8 -o solvated
package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.18 -cation POT -o ionized
quit
