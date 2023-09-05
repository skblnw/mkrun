mol new mol.pdb

set sel [atomselect top "resname ILE and name CD"]
$sel set name CD1
set sel [atomselect top "resname ASN and name OC1"]
$sel set name O
set sel [atomselect top "resname ASN and name OC2"]
$sel set name OXT
set sel [atomselect top "not hydrogen"]
$sel writepdb mol.pdb
quit
