package require mdff
#mdff griddx -i ../6tetramer_clip_rotated.mrc -o 6tetramer-grid.dx
#mdff gridpdb -psf ionized.psf -pdb ionized.pdb -o ionized-grid.pdb
#package require ssrestraints 
#ssrestraints -psf ionized.psf -pdb ionized.pdb -o ionized-extrabonds.txt -hbonds
#mol new ionized.psf	 
#mol addfile ionized.pdb
#package require cispeptide
#cispeptide restrain -o ionized-extrabonds-cispeptide.txt	 
#package require chirality
#chirality restrain -o ionized-extrabonds-chirality.txt

mdff setup -pbc -o 6tetramer -psf ionized.psf -pdb ionized.pdb -griddx 6tetramer-grid.dx -gridpdb ionized-grid.pdb -extrab {ionized-extrabonds.txt ionized-extrabonds-cispeptide.txt ionized-extrabonds-chirality.txt} -gscale 0.1 -numsteps 100000
