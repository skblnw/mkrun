#!/usr/bin/tclsh
#################################
## TCL script for writepdb from selection
## ukevi Apr 2014
## USAGE: vmd -dispdev text -e writepdb.tcl
#################################

mol new ../ionized.pdb type pdb waitfor all
set all [atomselect top "all"]

$all set beta 0
set sel [atomselect top "backbone or name P"]
$sel set beta 1
$all writepdb cons_bb_and_P.pdb

$all set beta 0
set sel [atomselect top "name CA P"]
$sel set beta 1
$all writepdb cons_CA_and_P.pdb

$all set beta 0
$all set occupancy 0
set sel [atomselect top "segname S1 S2 O1"]
$sel set beta 1
set sel [atomselect top "protein"]
$sel set occupancy 1
$all writepdb fix_memb.pdb

$all set beta 0
set sel [atomselect top "not {water or ion}"]
$sel set beta 1
$all writepdb fix_solute.pdb

$all set beta 0
$all set occupancy 0
set sel [atomselect top "protein"]
$sel set beta 1
$all writepdb ref_protein.pdb

quit
