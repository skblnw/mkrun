# If this script is used for published simulations, 
# please cite Balusek et al. Accelerating membrane 
# simulations with Hydrogen Mass Repartitioning. 
# J. Chem. Theory Comput. 2019.

# for background, see the paper
# Hopkins et al. Long-Time-Step Molecular Dynamics
# through Hydrogen Mass Repartitioning.
# J. Chem. Theory Comput. 11:1864-1874 (2015)
# https://pubs.acs.org/doi/full/10.1021/ct5010406

# checked against code from parmed here:
# https://parmed.github.io/ParmEd/html/_modules/parmed/tools/actions.html#HMassRepartition 

# forces VMD to respect atom names > 4 characters
set env(VMDNOMANGLEATOMNAMES) 1

set HMASS 3.024
set dowat 0

if { [llength $argv] != 2 } {
  puts stderr "\n This script repartitions the mass between heavy atoms and hydrogens in order
  to run simulations using a 4-fs time step. To run, type \n 

  vmd -dispdev text -e hmr.tcl -args <PSF> <PDB> \n"
  exit
}

set psf [lindex $argv 0]
set pdb [lindex $argv 1] 

mol new $psf
mol addfile $pdb

# don't need to do waters as SETTLE algorithm
# is stable without; SHAKE is not, however.
# note: if tested in NAMD, "useSettle on" is
# the default for waters with rigidBonds. 

if { $dowat } {
  set sel [atomselect top "resname TIP3 and name H1 H2"]
  $sel set mass $HMASS
  $sel delete 

  set sel [atomselect top "resname TIP3 and name OH2"]
  $sel set mass [expr 15.9994 - 2*($HMASS - 1.008)]
  $sel delete
}

# do the solute; there might be a faster
# way, but this seems most reliable.

# assume all atoms of mass 1 are hydrogens
set sel [atomselect top "(mass < 1.1 and mass > 0.99) and not water"]
foreach ind [$sel get index] {

  set temp [atomselect top "index $ind"]
  $temp set mass $HMASS
  set bondList [$temp getbonds]
  if { [llength $bondList] != 1 } {
    puts stderr "Hydrogen $ind has [llength $bondList] bond(s)!"
    exit
  }
  set temp2 [atomselect top "index $bondList"]
  $temp2 set mass [expr [$temp2 get mass] - ($HMASS - 1.008)]

  $temp delete
  $temp2 delete
}
$sel delete

# check masses
set sel [atomselect top "mass < 1"]
if { [$sel num] > 0 } {
  puts stderr "The mass of atoms [$sel get index] is less than 1!!!"
  exit
}
$sel delete 

set sel [atomselect top all]
$sel writepsf [string range $psf 0 end-4].hmr.psf 

quit

